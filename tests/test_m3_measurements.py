#!/usr/bin/env python3
"""The three M3-0 measurements (scripts/m3_0a, m3_0b): the k-mer encoding and the
k-mer width of a model. They decide the design of the unitig/k-mer article — whether a
unitig's k-mers carry its presence vector, and how wide the k-mer arm is — so the
arithmetic needs a fixture with a hand-computed answer."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
pytestmark = pytest.mark.unit


def _load(name: str):
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location(f"amrtest_{name}",
                                                  PROJECT_ROOT / "scripts" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


A = _load("m3_0a_unitig_kmer_homogeneity")
B = _load("m3_0b_kmer_feature_count")


def _codes(seq: str, k: int) -> np.ndarray:
    return A.canonical_codes(A._LUT[np.frombuffer(seq.encode(), dtype=np.uint8)], k)


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def test_a_sequence_yields_one_code_per_position():
    assert _codes("ACGTACGT", 4).size == 5
    assert _codes("ACGT", 4).size == 1
    assert _codes("ACG", 4).size == 0


def test_the_code_is_canonical_so_both_strands_agree():
    """A genome and its reverse strand must give the same k-mer, or a unitig called on
    one strand would look absent from a genome assembled on the other."""
    seq = "ACGTTGCAAGGCT"
    assert np.array_equal(_codes(seq, 5), _codes(_revcomp(seq), 5)[::-1])


def test_a_palindrome_and_its_complement_share_one_code():
    """ACGT is its own reverse complement, so forward and reverse pack to one value."""
    assert _codes("ACGT", 4)[0] == _codes("ACGT", 4)[0]
    assert _codes("AACGTT", 6)[0] == _codes(_revcomp("AACGTT"), 6)[0]
    assert _codes("AAAA", 4)[0] == _codes("TTTT", 4)[0]


def test_distinct_kmers_get_distinct_codes():
    codes = _codes("ACGTA", 4)          # ACGT and CGTA
    assert codes.size == 2 and codes[0] != codes[1]


def test_a_run_breaks_on_anything_that_is_not_acgt():
    """An N must not produce a k-mer that spans it, and a run shorter than k is dropped."""
    v = A._LUT[np.frombuffer(b"ACGTNACG", dtype=np.uint8)]
    runs = [r.size for r in A._acgt_runs(v, 4)]
    assert runs == [4]                  # ACGT kept, ACG too short
    v2 = A._LUT[np.frombuffer(b"acgtNacgtacgt", dtype=np.uint8)]
    assert [r.size for r in A._acgt_runs(v2, 4)] == [4, 8]


def test_unitig_lengths_counts_every_line_including_across_blocks(tmp_path):
    """The store holds tens of millions of lines, so the lengths come from newline
    positions per block; a small block forces the carry path."""
    (tmp_path / "sequences.txt").write_bytes(b"ACGT\nAC\nACGTACGTAC\n")
    for block in (1 << 24, 4, 3, 1):
        got = B.unitig_lengths(tmp_path, block=block)
        assert got.tolist() == [4, 2, 10], f"block={block}"


def test_unitig_lengths_keeps_a_last_line_without_a_newline(tmp_path):
    (tmp_path / "sequences.txt").write_bytes(b"ACGT\nACGTA")
    assert B.unitig_lengths(tmp_path, block=3).tolist() == [4, 5]


def test_the_kmer_width_of_a_pair_is_the_sum_over_kept_unitigs(tmp_path):
    """Three unitigs of length 31, 40 and 30 at k = 31 hold 1, 10 and 0 k-mers. The
    unitig arm counts one column per pattern; the k-mer arm counts one per k-mer, and
    nnz follows the carriers of the pattern each unitig belongs to."""
    lengths = np.asarray([31, 40, 30], dtype=np.int64)
    (tmp_path / "matrix_summary.json").write_text(
        '{"n_genomes": 100, "min_support": 10, "n_patterns": 2}')
    np.save(tmp_path / "members.npy", np.asarray([[0, 0], [1, 0], [2, 1]]))
    pd.DataFrame({"pattern_id": [0, 1], "n_present": [30, 20]}).to_csv(
        tmp_path / "patterns.csv", index=False)

    got = B.pair_counts(tmp_path, lengths, k=31)
    assert got["n_unitigs_kept"] == 3
    assert got["n_kmer_features"] == 1 + 10 + 0
    assert got["nnz_kmer"] == 1 * 30 + 10 * 30 + 0 * 20
    assert got["nnz_unitig"] == 30 + 20
    assert got["n_unitigs_below_k"] == 1
    assert got["median_unitig_length"] == 31.0


def test_a_members_file_from_another_run_is_refused(tmp_path):
    """members.npy pointing past the store means the matrix and the store are from
    different runs; a silent clip would report a k-mer width that is simply wrong."""
    (tmp_path / "matrix_summary.json").write_text(
        '{"n_genomes": 10, "min_support": 10, "n_patterns": 1}')
    np.save(tmp_path / "members.npy", np.asarray([[5, 0]]))
    pd.DataFrame({"pattern_id": [0], "n_present": [3]}).to_csv(
        tmp_path / "patterns.csv", index=False)
    with pytest.raises(IndexError, match="different runs"):
        B.pair_counts(tmp_path, np.asarray([31, 31], dtype=np.int64), k=31)

#!/usr/bin/env python3
"""Canonical k-mers and the organism k-mer store (lib.kmers) against a plain-Python count."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import kmers, matrix_store  # noqa: E402

pytestmark = pytest.mark.unit

_COMP = str.maketrans("ACGT", "TGCA")


def canonical(s: str) -> str:
    """Canonical form of a k-mer: the smaller of it and its reverse complement (A<C<G<T,
    which is also the order of the 2-bit codes)."""
    r = s.translate(_COMP)[::-1]
    return min(s, r)


def plain_kmers(seq: str, k: int) -> list[str]:
    return [canonical(seq[i:i + k]) for i in range(len(seq) - k + 1)
            if set(seq[i:i + k]) <= set("ACGT")]


def test_codes_match_a_plain_count_and_decode_back():
    rng = np.random.default_rng(1)
    seq = "".join(rng.choice(list("ACGT"), 300))
    seq = seq[:40] + "N" + seq[41:150] + "acgt" + seq[154:]           # a break, soft-masking
    k = 7
    got = kmers.codes(seq.encode(), k)
    want = plain_kmers(seq.upper(), k)
    assert got.size == len(want)
    assert [bytes(r).decode() for r in kmers.decode(got, k)] == want
    # a k-mer and its reverse complement give one code
    s = "ACCGTTA"
    assert kmers.codes(s.encode(), 7)[0] == kmers.codes(s.translate(_COMP)[::-1].encode(), 7)[0]


def test_no_kmer_spans_a_record_join_and_short_input_gives_none():
    assert kmers.codes(b"ACG", 5).size == 0
    assert kmers.codes(b"ACGTNACGT", 5).size == 0
    assert kmers.codes(b"ACGTA", 5).size == 1


@pytest.fixture
def organism(tmp_path):
    """12 assemblies of two records each, a unitig store whose presence table holds
    their sequences (so every k-mer is in a unitig), k = 9."""
    rng = np.random.default_rng(3)
    blocks = ["".join(rng.choice(list("ACGT"), 30)) for _ in range(14)]
    genomes = [f"g{i}" for i in range(12)]
    gdir = tmp_path / "genomes"
    gdir.mkdir()
    seqs = {}
    for i, g in enumerate(genomes):
        chosen = [b for j, b in enumerate(blocks) if j < 4 or rng.random() < 0.5]
        rec1 = "".join(chosen[: len(chosen) // 2])
        rec2 = "".join(chosen[len(chosen) // 2:]).lower() if i % 3 == 0 else \
            "".join(chosen[len(chosen) // 2:])
        seqs[g] = (rec1, rec2)
        (gdir / f"{g}.fna").write_text(f">{g}_1\n{rec1[:25]}\n{rec1[25:]}\n>{g}_2\n{rec2}\n")
    store = tmp_path / "unitig_store"
    (store / "call").mkdir(parents=True)
    pd.DataFrame({"Genome ID": genomes}).to_csv(store / "genomes.csv", index=False)
    (store / "store_summary.json").write_text('{"source_sha256": "abc"}')
    with open(store / "call" / "unitigs.rtab", "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(genomes) + "\n")
        for g in genomes:
            for rec in seqs[g]:
                f.write(rec.upper() + "\t" + "\t".join("1" for _ in genomes) + "\n")
    return store, [gdir / f"{g}.fna" for g in genomes], seqs


@pytest.mark.parametrize("workers", [1, 2])
def test_store_matches_a_plain_count(organism, tmp_path, workers):
    store_dir, paths, seqs = organism
    k, ms = 9, 2
    s = kmers.build_store(store_dir, paths, tmp_path / f"kmer_store{workers}", k=k,
                          min_support=ms, workers=workers, block_bytes=3)
    carriers: dict[str, set[int]] = {}
    for gi, g in enumerate(seqs):
        for rec in seqs[g]:
            for km in plain_kmers(rec.upper(), k):
                carriers.setdefault(km, set()).add(gi)
    n = len(seqs)
    kept = sorted(km for km, c in carriers.items() if ms <= len(c) <= n - ms)
    st = matrix_store.Store(tmp_path / f"kmer_store{workers}")
    assert st.genomes == list(seqs)
    assert s["n_kmers_graph"] == len(carriers) and s["n_kmers_in_no_genome"] == 0
    assert st.n_unitigs == len(kept)
    got = (tmp_path / f"kmer_store{workers}" / "sequences.txt").read_text().split()
    assert got == kept                  # codes ascending = strings ascending (A<C<G<T)
    bits = np.unpackbits(np.asarray(st.presence), axis=1)[:, :n]
    for row, km in zip(bits, got, strict=True):
        assert set(np.flatnonzero(row).tolist()) == carriers[km]
    assert not (tmp_path / f"kmer_store{workers}" / "presence_by_genome.bin").exists()


def test_store_refuses_missing_assemblies(organism, tmp_path):
    store_dir, paths, _ = organism
    with pytest.raises(FileNotFoundError):
        kmers.build_store(store_dir, paths[:-1] + [tmp_path / "nope.fna"], tmp_path / "x",
                          k=9, min_support=2)
    with pytest.raises(ValueError):
        kmers.build_store(store_dir, paths[:-1], tmp_path / "x", k=9, min_support=2)

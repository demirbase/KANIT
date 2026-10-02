#!/usr/bin/env python3
"""Organism store and model matrix (lib.matrix_store) on a small synthetic Rtab."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import matrix_store as ms  # noqa: E402

pytestmark = pytest.mark.unit

GENOMES = [f"g{i}" for i in range(12)]


def _write_rtab(path, rows):
    with open(path, "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(GENOMES) + "\n")
        for seq, bits in rows:
            f.write(seq + "\t" + "\t".join(str(b) for b in bits) + "\n")


@pytest.fixture
def rows():
    rng = np.random.default_rng(0)
    base = [list(rng.integers(0, 2, 12)) for _ in range(6)]
    base[0] = [1, 0] * 6                              # passes the filter in the subset below
    out = [(f"U{i:02d}" + "ACGT" * (i % 3), b) for i, b in enumerate(base)]
    out.append(("TWIN0", list(base[0])))             # identical to U00 -> same pattern
    out.append(("COMPL0", [1 - b for b in base[0]]))  # complement -> its own pattern
    out.append(("RARE", [1] + [0] * 11))              # present in 1 genome
    out.append(("CORE", [1] * 12))                    # present in all genomes
    out.append(("NEARCORE", [0] + [1] * 11))          # absent from 1 genome
    return out


def test_build_store_keeps_rows_that_can_pass_a_subset(tmp_path, rows):
    rtab = tmp_path / "u.rtab"
    _write_rtab(rtab, rows)
    s = ms.build_store(rtab, tmp_path / "store", min_support=2, block_bytes=64)  # tiny blocks
    store = ms.Store(tmp_path / "store")
    kept = [(q, b) for q, b in rows if 2 <= sum(b) <= 10]
    assert s["n_unitigs_seen"] == len(rows) and store.n_unitigs == len(kept)
    assert store.genomes == GENOMES
    got = np.unpackbits(np.asarray(store.presence), axis=1)[:, :12]
    assert got.tolist() == [list(b) for _, b in kept]
    assert store.sequences(range(store.n_unitigs)) == [q for q, _ in kept]


def test_build_store_rejects_a_malformed_line(tmp_path):
    rtab = tmp_path / "bad.rtab"
    with open(rtab, "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(GENOMES) + "\n")
        f.write("AAA\t" + "\t".join(["1"] * 11) + "\n")   # one value short
    with pytest.raises(ValueError, match="malformed"):
        ms.build_store(rtab, tmp_path / "store", min_support=1)


def test_model_matrix_filters_and_collapses_patterns(tmp_path, rows):
    rtab = tmp_path / "u.rtab"
    _write_rtab(rtab, rows)
    ms.build_store(rtab, tmp_path / "store", min_support=1)
    sub = [9, 1, 2, 3, 4, 5, 6, 7, 8, 10]              # 10 genomes, not in store order
    genomes = pd.DataFrame({"Genome ID": [GENOMES[i] for i in sub],
                            "label": [i % 2 for i in sub], "lineage": [i // 4 for i in sub]})
    s = ms.build_model_matrix(tmp_path / "store", genomes, tmp_path / "model",
                              min_support=2, block_unitigs=3)          # tiny blocks
    m = ms.ModelMatrix(tmp_path / "model")
    store = ms.Store(tmp_path / "store")
    pres = np.unpackbits(np.asarray(store.presence), axis=1)[:, sub]   # (unitigs, 10)

    # symmetric frequency filter over the model's genomes
    expect_kept = [i for i in range(store.n_unitigs) if 2 <= pres[i].sum() <= 8]
    mem = m.members()
    assert sorted(mem["unitig_index"]) == expect_kept
    # identical columns share a pattern, a complement does not
    seqs = store.sequences(range(store.n_unitigs))
    pid = dict(zip(mem["unitig_index"], mem["pattern_id"], strict=True))
    i_u00, i_twin, i_comp = seqs.index("U00"), seqs.index("TWIN0"), seqs.index("COMPL0")
    assert pid[i_u00] == pid[i_twin] and pid[i_u00] != pid[i_comp]
    # X holds each pattern once, rows in the order given
    dense = m.rows(range(m.n_genomes))
    assert dense.shape == (10, s["n_patterns"])
    for u, p in pid.items():
        assert dense[:, p].tolist() == pres[u].tolist()
        assert m.pattern(p).tolist() == pres[u].tolist()
    assert len({tuple(dense[:, p]) for p in range(m.n_patterns)}) == m.n_patterns
    # batches return the requested rows
    pick = [7, 0, 3]
    assert np.vstack(list(m.row_batches(pick, 2))).tolist() == dense[pick].tolist()
    assert m.labels.tolist() == genomes["label"].tolist()


def test_model_matrix_needs_every_genome_in_the_store(tmp_path, rows):
    rtab = tmp_path / "u.rtab"
    _write_rtab(rtab, rows)
    ms.build_store(rtab, tmp_path / "store", min_support=1)
    genomes = pd.DataFrame({"Genome ID": ["g1", "gX"], "label": [0, 1], "lineage": [1, 1]})
    with pytest.raises(ValueError, match="not in the store"):
        ms.build_model_matrix(tmp_path / "store", genomes, tmp_path / "model", min_support=1)

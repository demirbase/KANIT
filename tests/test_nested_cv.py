#!/usr/bin/env python3
"""Training layer (lib.train) and the nested-CV runner (04_nested_cv.py) on a
small synthetic model matrix with a planted resistance signal."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

pytest.importorskip("xgboost")
pytest.importorskip("optuna")

from lib import folds, matrix_store, train  # noqa: E402
from lib.lineage import no_group_leakage  # noqa: E402

pytestmark = pytest.mark.unit

CV = {"n_repeats": 1, "n_folds": 3, "min_minority_per_test_fold": 5, "max_seed_attempts": 20,
      "threshold": 0.5, "n_bootstrap": 30, "reliability_bins": 5}
HPO = {"n_trials": 3, "search_max_genomes": 80, "pruner_startup_trials": 1,
       "pruner_warmup_rounds": 5, "early_stopping_rounds": 10, "max_rounds": 100,
       "batch_rows": 50}


@pytest.fixture(scope="module")
def model(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("cv")
    rng = np.random.default_rng(0)
    n = 160
    ids = [f"g{i:03d}" for i in range(n)]
    lineage = np.arange(n) // 5                       # 32 lineages of 5
    signal = rng.random(n) < 0.4
    y = np.where(rng.random(n) < 0.9, signal, ~signal).astype(int)
    rows = [("SIGNAL", signal.astype(int))]
    rows += [(f"U{j}", (rng.random(n) < rng.uniform(0.3, 0.7)).astype(int)) for j in range(200)]
    rtab = tmp / "u.rtab"
    with open(rtab, "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(ids) + "\n")
        for seq, bits in rows:
            f.write(seq + "\t" + "\t".join(map(str, bits)) + "\n")
    matrix_store.build_store(rtab, tmp / "store", min_support=5)
    genomes = pd.DataFrame({"Genome ID": ids, "label": y, "lineage": lineage})
    matrix_store.build_model_matrix(tmp / "store", genomes, tmp / "model", min_support=5)
    return matrix_store.ModelMatrix(tmp / "model"), tmp


def test_search_space_and_weights():
    assert train.colsample_range(1_000_000) == pytest.approx((5e-4, 0.02))
    assert train.colsample_range(100) == pytest.approx((0.05, 1.0))
    w = train.class_weight([1, 0, 0, 0])
    assert w.tolist() == [3.0, 1.0, 1.0, 1.0]
    with pytest.raises(ValueError):
        train.class_weight([0, 0])


def test_search_subsample_keeps_the_class_ratio_and_is_seeded():
    rng = np.random.default_rng(1)
    y = (rng.random(3000) < 0.25).astype(int)
    groups = rng.integers(0, 50, 3000)
    a = train.search_subsample(y, groups, 1000, seed=7)
    assert len(a) == 1000 and len(set(a.tolist())) == 1000
    assert abs(y[a].mean() - y.mean()) < 0.01
    assert (train.search_subsample(y, groups, 1000, seed=7) == a).all()
    assert len(train.search_subsample(y[:500], groups[:500], 1000, seed=7)) == 500


def test_inner_split_matches_the_arm(model):
    mm, _ = model
    y, groups = mm.labels, mm.genomes["lineage"].to_numpy()
    tr, va = train.inner_split(y, groups, folds.LINEAGE_AWARE, seed=3)
    m_tr = np.zeros(len(y), bool)
    m_va = np.zeros(len(y), bool)
    m_tr[tr], m_va[va] = True, True
    assert no_group_leakage(m_tr, m_va, groups)
    assert set(y[va].tolist()) == {0, 1} and set(y[tr].tolist()) == {0, 1}


def test_runner_end_to_end(model, load_script):
    mm, tmp = model
    r = load_script("04_nested_cv.py")
    out = tmp / "out"
    design = r.run_folds(mm, out, CV)
    assert design["evaluable"] and (out / "folds.csv").exists()
    seeds = pd.read_csv(out / "seeds.csv")
    units = [(s.arm, int(s.repeat), k) for s in seeds.itertuples() for k in range(3)]
    for u in units[:-1]:
        r.run_unit(mm, out, *u, HPO, threads=1)
    with pytest.raises(SystemExit, match="not finished"):        # a missing fold stops it
        r.run_metrics(mm, out, CV)
    r.run_unit(mm, out, *units[-1], HPO, threads=1)
    r.run_final(mm, out, HPO, threads=1)
    m = r.run_metrics(mm, out, CV)
    oof = pd.read_csv(out / "oof_predictions.csv")
    assert (oof.groupby(["arm", "repeat"])["genome_id"].nunique() == mm.n_genomes).all()
    for arm in folds.ARMS:
        assert m["arms"][arm]["roc_auc"] > 0.8               # the planted signal is found
    assert "lineage_blind_minus_lineage_aware" in m["bootstrap"]
    rec = json.loads((out / "units" / folds.unit_name(*units[0]) / "record.json").read_text())
    assert rec["n_trees"] >= 1 and rec["n_trials_complete"] + rec["n_trials_pruned"] == 3
    assert (out / "final" / "model.ubj").exists()

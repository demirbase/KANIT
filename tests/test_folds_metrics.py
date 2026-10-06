#!/usr/bin/env python3
"""Outer folds (lib.folds) and out-of-fold metrics (lib.oof_metrics)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import folds, oof_metrics  # noqa: E402
from lib.lineage import no_group_leakage  # noqa: E402

pytestmark = pytest.mark.unit


@pytest.fixture
def data():
    rng = np.random.default_rng(0)
    n = 400
    groups = rng.integers(0, 60, n)                     # 60 lineages
    y = (rng.random(n) < 0.3).astype(int)               # minority: resistant
    return y, groups


def test_both_arms_partition_and_the_aware_arm_keeps_lineages_whole(data):
    y, groups = data
    for arm in folds.ARMS:
        fold_of = folds.split(y, groups, arm, seed=1000, n_folds=5)
        assert sorted(set(fold_of.tolist())) == [0, 1, 2, 3, 4]
        if arm == folds.LINEAGE_AWARE:
            for k in range(5):
                assert no_group_leakage(fold_of != k, fold_of == k, groups)


def _dominated():
    """One lineage holds 40% of the genomes, as in the pilot (A. baumannii)."""
    rng = np.random.default_rng(1)
    sizes = np.array([389, 60, 45, 40, 30, 25, 20, 15, 12, 10] + [2] * 20 + [1] * 300)
    groups = np.repeat(np.arange(len(sizes)), sizes)
    p_r = np.r_[0.65, rng.uniform(0.05, 0.9, len(sizes) - 1)]
    return (rng.random(len(groups)) < p_r[groups]).astype(int), groups


def test_the_repeats_of_the_lineage_aware_arm_differ():
    """The seed changes the folds of the lineages besides the largest, not only the order
    of ties (with scikit-learn's StratifiedGroupKFold the repeats had nearly the same
    folds)."""
    from sklearn.metrics import adjusted_rand_score

    y, groups = _dominated()
    rest = groups != 0
    f1, f2 = (folds.split(y, groups, folds.LINEAGE_AWARE, s, 5) for s in (1000, 2000))
    assert (folds.split(y, groups, folds.LINEAGE_AWARE, 1000, 5) == f1).all()   # reproducible
    assert adjusted_rand_score(f1[rest], f2[rest]) < 0.2
    for f in (f1, f2):
        for k in range(5):
            assert no_group_leakage(f != k, f == k, groups)
        assert (f == f[groups == 0][0]).sum() == (groups == 0).sum()   # the largest alone


def test_the_inner_validation_fold_is_about_a_fifth():
    from lib import train

    y, groups = _dominated()
    _, va = train.inner_split(y, groups, folds.LINEAGE_AWARE, seed=0)
    assert not (groups[va] == 0).any()                 # not the largest lineage's fold
    assert 0.1 < len(va) / len(y) < 0.3


def test_seed_sequence_and_the_class_balance_rule(data):
    y, groups = data
    ok = folds.assign_repeat(y, groups, folds.LINEAGE_AWARE, 3, n_folds=5,
                             min_minority=1, max_attempts=20)
    assert ok.seed == 3000 and ok.attempts == 1 and ok.evaluable
    # a rule no split can meet: every one of the 20 seeds is tried, none is kept
    bad = folds.assign_repeat(y, groups, folds.LINEAGE_AWARE, 3, n_folds=5,
                              min_minority=10_000, max_attempts=20)
    assert not bad.evaluable and bad.seed is None and bad.attempts == 20


def test_fewer_lineages_than_folds_is_not_evaluable():
    y = np.array([0, 1] * 10)
    groups = np.array([0, 1, 2] * 6 + [0, 1])
    s = folds.assign_repeat(y, groups, folds.LINEAGE_AWARE, 1, n_folds=5,
                            min_minority=1, max_attempts=20)
    assert not s.evaluable


def test_tables_record_every_split(data):
    y, groups = data
    splits = folds.assign_all(y, groups, n_repeats=2, n_folds=5, min_minority=5, max_attempts=20)
    ids = [f"g{i}" for i in range(len(y))]
    ft = folds.folds_table(ids, splits)
    assert len(ft) == 2 * 2 * len(y)                     # 2 arms x 2 repeats x every genome
    comp = folds.composition_table(y, groups, splits, 5)
    assert (comp.groupby(["arm", "repeat"])["n"].sum() == len(y)).all()
    assert (comp["n_resistant"] + comp["n_susceptible"] == comp["n"]).all()
    assert set(folds.seeds_table(splits)["seed"]) == {1000, 2000}


def _oof(y, p, repeats=2, n_folds=2):
    n = len(y)
    return pd.concat([pd.DataFrame({"repeat": r, "fold": np.arange(n) % n_folds,
                                    "genome_id": [f"g{i}" for i in range(n)], "y": y, "p": p})
                      for r in range(1, repeats + 1)], ignore_index=True)


def test_pooled_metrics_and_a_single_class_fold():
    y = np.array([0, 1, 0, 1, 0, 1, 0, 0])
    oof = _oof(y, y * 0.9 + 0.05)
    s = oof_metrics.summary(oof)
    assert s["roc_auc"] == 1.0 and s["balanced_accuracy"] == 1.0 and s["n_repeats"] == 2
    # fold 0 holds genomes 0, 2, 4, 6: all susceptible -> undefined AUC, still pooled
    pf = oof_metrics.per_fold(oof)
    assert pf.loc[pf["fold"] == 0, "roc_auc"].isna().all()
    assert s["n_folds_single_class"] == 2


def test_cluster_bootstrap_is_paired_across_arms():
    rng = np.random.default_rng(1)
    n = 300
    y = (rng.random(n) < 0.4).astype(int)
    p = np.clip(y * 0.6 + rng.random(n) * 0.5, 0, 1)
    oof = _oof(y, p)
    lineage_of = {f"g{i}": i % 40 for i in range(n)}
    ci = oof_metrics.cluster_bootstrap({"a": oof, "b": oof.copy()}, lineage_of,
                                       n_boot=200, seed=0, diff=("a", "b"))
    point = oof_metrics.summary(oof)["roc_auc"]
    assert ci["a"]["low"] <= point <= ci["a"]["high"]
    # identical arms under the same resamples: the difference is exactly zero
    d = ci["a_minus_b"]
    assert d["low"] == 0.0 and d["high"] == 0.0 and ci["n_lineages"] == 40


def test_reliability_bins_cover_all_predictions():
    y = np.array([0, 1, 1, 0, 1])
    p = np.array([0.05, 0.95, 1.0, 0.0, 0.5])
    rel = oof_metrics.reliability(_oof(y, p, repeats=1), n_bins=10)
    assert rel["n"].sum() == 5 and len(rel) == 10

#!/usr/bin/env python3
"""Label permutation test (lib.label_permutation, 12b_label_permutation.py) on a
small synthetic model with real outer-fold models of 04."""
import importlib.util
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

from lib import folds, matrix_store  # noqa: E402
from lib import label_permutation as lp  # noqa: E402

pytestmark = pytest.mark.unit

CV = {"n_repeats": 1, "n_folds": 3, "min_minority_per_test_fold": 5, "max_seed_attempts": 20,
      "threshold": 0.5, "n_bootstrap": 30, "reliability_bins": 5}
HPO = {"n_trials": 3, "search_max_genomes": 80, "pruner_startup_trials": 1,
       "pruner_warmup_rounds": 5, "early_stopping_rounds": 10, "max_rounds": 100,
       "batch_rows": 50}
CFG = {"n_permutations": 12, "seed": 0, "chunk": 5, "alpha": 0.05}


def _script(name):
    spec = importlib.util.spec_from_file_location("amrtest_" + Path(name).stem,
                                                  PROJECT_ROOT / "scripts" / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def cv(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("lp")
    rng = np.random.default_rng(0)
    n = 150
    ids = [f"g{i:03d}" for i in range(n)]
    signal = rng.random(n) < 0.4
    y = np.where(rng.random(n) < 0.9, signal, ~signal).astype(int)
    rows = [("SIGNAL", signal.astype(int))]
    rows += [(f"U{j}", (rng.random(n) < rng.uniform(0.3, 0.7)).astype(int)) for j in range(40)]
    with open(tmp / "u.rtab", "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(ids) + "\n")
        for seq, bits in rows:
            f.write(seq + "\t" + "\t".join(map(str, bits)) + "\n")
    matrix_store.build_store(tmp / "u.rtab", tmp / "store", min_support=5)
    genomes = pd.DataFrame({"Genome ID": ids, "label": y, "lineage": np.arange(n) // 5})
    matrix_store.build_model_matrix(tmp / "store", genomes, tmp / "model", min_support=5)
    mm = matrix_store.ModelMatrix(tmp / "model")
    r = _script("04_nested_cv.py")
    assert r.run_folds(mm, tmp / "cv", CV)["evaluable"]
    for k in range(CV["n_folds"]):
        r.run_unit(mm, tmp / "cv", folds.LINEAGE_AWARE, 1, k, HPO, threads=1)
    return mm, tmp


def test_permuted_labels_are_seeded_and_keep_the_counts():
    y = np.array([1] * 30 + [0] * 70)
    a, b = lp.permuted_labels(y, 0, 1), lp.permuted_labels(y, 0, 2)
    assert (a == lp.permuted_labels(y, 0, 1)).all() and not (a == b).all()
    assert a.sum() == 30 and b.sum() == 30
    assert [c.tolist() for c in lp.chunks(12, 5)] == [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9], [10, 11]]


def test_refit_on_the_real_labels_reproduces_04(cv):
    """The refit uses the fold's own hyperparameters, trees, seed and class weights."""
    mm, tmp = cv
    fold_list = lp.folds_of(mm, tmp / "cv")
    assert [f.fold for f in fold_list] == [0, 1, 2]
    assert sorted(np.concatenate([f.test_rows for f in fold_list]).tolist()) == list(
        range(mm.n_genomes))
    for f in fold_list:
        refit = lp.Refitter(mm, f, threads=1)
        assert refit.predict(mm.labels.astype(int)) == pytest.approx(f.p, abs=1e-6)


def test_run_and_metrics(cv, tmp_path):
    mm, tmp = cv
    s12b = _script("12b_label_permutation.py")
    out = tmp_path / "lp"
    s12b.run(mm, tmp / "cv", out, CFG, fold=0, chunk=1, threads=1)
    assert [p.name for p in sorted((out / "null").iterdir())] == ["fold0_chunk001.npz"]
    with pytest.raises(FileNotFoundError, match="8 null chunk"):
        s12b.metrics(mm, tmp / "cv", out, CFG)
    before = (out / "null" / "fold0_chunk001.npz").stat().st_mtime_ns
    s12b.run(mm, tmp / "cv", out, CFG, threads=1)
    assert (out / "null" / "fold0_chunk001.npz").stat().st_mtime_ns == before   # not rerun
    s = s12b.metrics(mm, tmp / "cv", out, CFG)
    null = pd.read_csv(out / "label_permutation_null.csv")["auc"]
    assert len(null) == 12 and abs(null.mean() - 0.5) < 0.15
    assert s["auc_observed"] > 0.8 and s["p"] == pytest.approx(1 / 13) and s["z"] > 3
    assert json.loads((out / "label_permutation.json").read_text())["n_permutations"] == 12


def test_across_models_flags_the_non_significant():
    rows = pd.DataFrame({"organism": list("abcd"), "p": [0.0099, 0.0099, 0.5, np.nan]})
    t = lp.across_models(rows, 0.05).set_index("organism")
    assert t.loc["a", "q"] == pytest.approx(0.01485) and t.loc["a", "flag"] == ""
    assert t.loc["c", "flag"] == lp.FLAG
    assert np.isnan(t.loc["d", "q"]) and t.loc["d", "flag"] == ""


def test_across_reads_the_panel(tmp_path, monkeypatch):
    s12b = _script("12b_label_permutation.py")
    root = tmp_path / "r"
    config = {"paths_organism": {
        "panel_dir": str(root / "panel"), "cv_dir": str(root / "{organism}/{antibiotic}/cv"),
        "label_permutation_dir": str(root / "{organism}/{antibiotic}/lp"),
        "cross_model_dir": str(root / "cross")}, "label_permutation": CFG}
    (root / "panel").mkdir(parents=True)
    pd.DataFrame({"organism": ["ec", "ec", "kp", "kp"], "antibiotic": ["amp", "cip", "amp", "x"],
                  "decision": ["included", "included", "included", "excluded"]}).to_csv(
        root / "panel" / "panel_decisions.csv", index=False)
    for org, ab, evaluable in (("ec", "amp", True), ("ec", "cip", True), ("kp", "amp", False)):
        (root / org / ab / "cv").mkdir(parents=True)
        (root / org / ab / "cv" / "cv_design.json").write_text(json.dumps({"evaluable": evaluable}))
    summary = {"auc_observed": 0.9, "null_mean": 0.5, "null_sd": 0.03, "z": 13.3,
               "n_permutations": 100, "p": 1 / 101}
    (root / "ec" / "amp" / "lp").mkdir(parents=True)
    (root / "ec" / "amp" / "lp" / "label_permutation.json").write_text(json.dumps(summary))
    with pytest.raises(SystemExit, match="ec/cip"):
        s12b.across(config)
    (root / "ec" / "cip" / "lp").mkdir(parents=True)
    (root / "ec" / "cip" / "lp" / "label_permutation.json").write_text(
        json.dumps({**summary, "auc_observed": 0.52, "z": 0.4, "p": 0.3}))
    t = s12b.across(config).set_index(["organism", "antibiotic"])
    assert (root / "cross" / "label_permutation.csv").exists()
    assert t["status"].to_dict() == {("ec", "amp"): "tested", ("ec", "cip"): "tested",
                                     ("kp", "amp"): "not_evaluable"}
    assert t["flag"].to_dict() == {("ec", "amp"): "", ("ec", "cip"): lp.FLAG, ("kp", "amp"): ""}

#!/usr/bin/env python3
"""Prevalence (lib.prevalence, 10_prevalence.py) and MDA (lib.mda, 12_mda.py) on
small synthetic model matrices; the MDA uses real outer-fold models of 04."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control, fisher_exact
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

pytest.importorskip("xgboost")
pytest.importorskip("optuna")

from lib import folds, matrix_store, mda  # noqa: E402
from lib.prevalence import prevalence_layer  # noqa: E402

pytestmark = pytest.mark.unit

CV = {"n_repeats": 2, "n_folds": 3, "min_minority_per_test_fold": 5, "max_seed_attempts": 20,
      "threshold": 0.5, "n_bootstrap": 30, "reliability_bins": 5}
HPO = {"n_trials": 3, "search_max_genomes": 80, "pruner_startup_trials": 1,
       "pruner_warmup_rounds": 5, "early_stopping_rounds": 10, "max_rounds": 100,
       "batch_rows": 50}
R = 200


def _model(tmp, rows, labels, lineage, min_support):
    ids = [f"g{i:03d}" for i in range(len(labels))]
    rtab = tmp / "u.rtab"
    with open(rtab, "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(ids) + "\n")
        for seq, bits in rows:
            f.write(seq + "\t" + "\t".join(map(str, bits)) + "\n")
    matrix_store.build_store(rtab, tmp / "store", min_support=min_support)
    genomes = pd.DataFrame({"Genome ID": ids, "label": labels, "lineage": lineage})
    matrix_store.build_model_matrix(tmp / "store", genomes, tmp / "model",
                                    min_support=min_support)
    mm = matrix_store.ModelMatrix(tmp / "model")
    store = matrix_store.Store(tmp / "store")
    m = mm.members()
    seq = dict(zip(m["unitig_index"], store.sequences(m["unitig_index"]), strict=True))
    return mm, {seq[u]: int(p) for u, p in zip(m["unitig_index"], m["pattern_id"], strict=True)}


def _script(name):
    spec = importlib.util.spec_from_file_location("amrtest_" + Path(name).stem,
                                                  PROJECT_ROOT / "scripts" / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def cv(tmp_path_factory):
    """A planted signal, a near copy of it (|r| > 0.9) and noise; 04's lineage-aware
    outer folds, 2 repeats × 3 folds."""
    tmp = tmp_path_factory.mktemp("mda")
    rng = np.random.default_rng(0)
    n = 160
    signal = rng.random(n) < 0.4
    y = np.where(rng.random(n) < 0.9, signal, ~signal).astype(int)
    near = signal.copy()
    near[rng.choice(n, 3, replace=False)] ^= True
    rows = [("SIGNAL", signal.astype(int)), ("NEAR", near.astype(int))]
    rows += [(f"U{j}", (rng.random(n) < rng.uniform(0.3, 0.7)).astype(int)) for j in range(60)]
    mm, pid = _model(tmp, rows, y, np.arange(n) // 5, 5)
    r = _script("04_nested_cv.py")
    out = tmp / "cv"
    assert r.run_folds(mm, out, CV)["evaluable"]
    for repeat in range(1, CV["n_repeats"] + 1):
        for fold in range(CV["n_folds"]):
            r.run_unit(mm, out, folds.LINEAGE_AWARE, repeat, fold, HPO, threads=1)
    design = mda.Design(mda.load_fold_models(mm, out), R, 0, mm.labels.astype(int))
    return mm, design, pid, tmp


# ---- helpers ---------------------------------------------------------------
def test_auc_rows_matches_sklearn_with_ties():
    rng = np.random.default_rng(1)
    y = rng.random(300) < 0.3
    s = np.round(rng.random((5, 300)), 1)                    # many ties
    expected = [roc_auc_score(y, row) for row in s]
    assert mda.auc_rows(s, y) == pytest.approx(expected, abs=1e-12)


def test_columns_read_several_patterns(cv):
    mm = cv[0]
    ids = [3, 0, 7]
    assert (mm.columns(ids) == np.stack([mm.pattern(i) for i in ids], axis=1)).all()
    with pytest.raises(IndexError):
        mm.columns([mm.n_patterns])


def test_correlation_clusters_link_chains_and_complements(tmp_path):
    a = np.array([1] * 20 + [0] * 20)
    b = a.copy()
    b[:1] = 0                                          # r(a, b) > 0.9
    c = b.copy()
    c[20:21] = 1                                       # r(b, c) > 0.9
    rows = [("A", a), ("B", b), ("C", c), ("NOT_A", 1 - a), ("FAR", np.tile([1, 0], 20))]
    mm, pid = _model(tmp_path, rows, np.tile([0, 1], 20), np.arange(40), 2)
    groups = mda.correlation_clusters(mm, sorted(pid.values()), 0.9)
    by = {name: next(i for i, g in enumerate(groups) if pid[name] in g) for name in pid}
    assert by["A"] == by["B"] == by["C"] == by["NOT_A"] != by["FAR"]


# ---- prevalence ------------------------------------------------------------
def test_prevalence_layer(tmp_path):
    labels = np.array([1] * 20 + [0] * 20)
    enriched = np.array([1] * 16 + [0] * 4 + [1] * 4 + [0] * 16)
    flat = np.tile([1, 0], 20)
    depleted = 1 - enriched
    rows = [("ENR", enriched), ("FLAT", flat), ("DEP", depleted)]
    mm, pid = _model(tmp_path, rows, labels, np.arange(40), 2)
    layer = prevalence_layer(mm, [pid["ENR"], pid["FLAT"], pid["DEP"]]).set_index("pattern_id")
    enr, fl, dep = layer.loc[pid["ENR"]], layer.loc[pid["FLAT"]], layer.loc[pid["DEP"]]
    assert (enr["present_resistant"], enr["present_susceptible"]) == (16, 4)
    assert enr["delta"] == pytest.approx(0.6) and enr["direction"] == "R" and enr["passes"]
    assert enr["fisher_p"] == pytest.approx(fisher_exact([[16, 4], [4, 16]])[1])
    assert fl["direction"] == "none" and not fl["passes"]
    assert dep["direction"] == "S" and dep["passes"]
    assert layer["q"].to_numpy() == pytest.approx(
        false_discovery_control(layer["fisher_p"].to_numpy(), method="bh"))


# ---- MDA -------------------------------------------------------------------
def test_observed_auc_is_the_pooled_auc_of_04(cv):
    _, design, _, tmp = cv
    for r in design.by_repeat:
        oof = pd.concat([pd.read_csv(tmp / "cv" / "units" /
                                     folds.unit_name(folds.LINEAGE_AWARE, r, k) / "oof.csv")
                         for k in range(CV["n_folds"])])
        assert design.auc[r] == pytest.approx(roc_auc_score(oof["y"], oof["p"]), abs=1e-12)


def test_permutation_lookup_matches_brute_force(cv):
    """Permuting the column and predicting again gives the same AUC in every draw."""
    mm, design, pid, _ = cv
    j = pid["SIGNAL"]
    fast, n_used = mda.permuted_auc(design, [j])
    assert n_used > 0
    for b in (0, 1, 57):
        aucs = []
        for r in sorted(design.by_repeat):
            ys, ps = [], []
            for i in design.by_repeat[r]:
                fm = design.models[i]
                x = mm.rows(fm.rows).astype(np.float32)
                x[:, j] = x[design.perms[i][b].astype(np.int64), j]
                ps.append(fm.booster.inplace_predict(x))
                ys.append(mm.labels[fm.rows])
            aucs.append(roc_auc_score(np.concatenate(ys), np.concatenate(ps)))
        assert fast[b] == pytest.approx(np.mean(aucs), abs=1e-12)


def test_mda_layer_finds_the_signal_and_its_cluster(cv):
    mm, design, pid, _ = cv
    unused = sorted(set(range(mm.n_patterns))
                    - set(np.concatenate([fm.used for fm in design.models]).tolist()))
    cands = [pid["SIGNAL"], pid["NEAR"], pid["U0"], pid["U1"], pid["U2"], *unused[:1]]
    main, clus = mda.mda_layer(mm, design, cands)
    by = main.set_index("pattern_id")
    sig = by.loc[pid["SIGNAL"]]
    # SIGNAL and its near copy share the importance (a fold model may split on either):
    # permuting SIGNAL alone shows part of it, permuting the cluster all of it
    assert sig["mda"] > 0.05 and sig["p"] == pytest.approx(1 / (R + 1)) and sig["passes"]
    for u in unused[:1]:                                   # never split on: no effect
        assert by.loc[u, "n_fold_models_using"] == 0
        assert by.loc[u, "mda"] == 0 and by.loc[u, "p"] == 1 and not by.loc[u, "passes"]
    c = clus.set_index("pattern_id")
    assert c.loc[pid["SIGNAL"], "cluster"] == c.loc[pid["NEAR"], "cluster"]
    assert c.loc[pid["SIGNAL"], "cluster_size"] == 2 and c.loc[pid["SIGNAL"], "passes"]
    assert c.loc[pid["SIGNAL"], "mda"] > 0.1
    again, _ = mda.mda_layer(mm, design, cands)
    pd.testing.assert_frame_equal(main, again)             # the draws are seeded
    none, n_used = mda.permuted_auc(design, [mm.n_patterns + 1])
    assert n_used == 0 and (none == design.observed).all()


def test_stored_model_must_reproduce_04(cv, tmp_path):
    mm, design, _, tmp = cv
    import shutil
    bad = tmp_path / "cv"
    shutil.copytree(tmp / "cv", bad)
    unit = bad / "units" / folds.unit_name(folds.LINEAGE_AWARE, 1, 0)
    oof = pd.read_csv(unit / "oof.csv")
    oof.assign(p=1 - oof["p"]).to_csv(unit / "oof.csv", index=False)
    with pytest.raises(ValueError, match="does not reproduce"):
        mda.load_fold_models(mm, bad)
    oof.to_csv(unit / "oof.csv", index=False)
    shutil.rmtree(bad / "units" / folds.unit_name(folds.LINEAGE_AWARE, 2, 1))
    with pytest.raises(FileNotFoundError):
        mda.load_fold_models(mm, bad)


def test_layer_scripts_end_to_end(cv, tmp_path, monkeypatch):
    mm, _, pid, tmp = cv
    cand = tmp_path / "candidates.csv"
    pd.DataFrame({"pattern_id": [pid["SIGNAL"], pid["U3"], pid["NEAR"]]}).to_csv(cand, index=False)
    config = {"paths_organism": {"candidates_file": str(cand), "matrix_dir": str(tmp / "model"),
                                 "cv_dir": str(tmp / "cv"), "layers_dir": str(tmp_path / "layers")},
              "prevalence": {"min_delta": 0.10, "alpha": 0.05},
              "mda": {"n_permutations": 50, "seed": 0, "alpha": 0.05, "cluster_r": 0.9}}
    monkeypatch.setattr(sys, "argv", ["x", "--organism", "o", "--antibiotic", "a"])
    for name in ("10_prevalence.py", "12_mda.py"):
        m = _script(name)
        monkeypatch.setattr(m, "load_config", lambda: config)
        m.main()
    layers = tmp_path / "layers"
    for f in ("prevalence.csv", "mda.csv"):
        d = pd.read_csv(layers / f)
        assert {"pattern_id", "passes"} <= set(d.columns) and len(d) == 3
    s = json.loads((layers / "mda_summary.json").read_text())
    assert s["n_fold_models"] == 6 and s["n_permutations"] == 50
    design_file = tmp_path / "cv_na" / "cv_design.json"
    design_file.parent.mkdir()
    design_file.write_text(json.dumps({"evaluable": False}))
    config["paths_organism"]["cv_dir"] = str(design_file.parent)
    with pytest.raises(SystemExit) as e:
        m.main()
    assert e.value.code == folds.NOT_EVALUABLE

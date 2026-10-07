#!/usr/bin/env python3
"""Stability selection (lib.cpss, 13_cpss.py) and the pyseer step (lib.pyseer_lmm,
14_pyseer.py) on a small synthetic model; CPSS uses the real final model of 04."""
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2, chi2_contingency

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

pytest.importorskip("xgboost")
pytest.importorskip("optuna")

from lib import cpss, matrix_store  # noqa: E402
from lib import pyseer_lmm as pl  # noqa: E402

pytestmark = pytest.mark.unit

CV = {"n_repeats": 1, "n_folds": 3, "min_minority_per_test_fold": 5, "max_seed_attempts": 20,
      "threshold": 0.5, "n_bootstrap": 30, "reliability_bins": 5}
HPO = {"n_trials": 3, "search_max_genomes": 80, "pruner_startup_trials": 1,
       "pruner_warmup_rounds": 5, "early_stopping_rounds": 10, "max_rounds": 100,
       "batch_rows": 50}
CPSS = {"prefilter": 30, "n_pairs": 6, "q": 5, "pi_threshold": 0.6, "seed": 0, "chunk": 4}


def _script(name):
    spec = importlib.util.spec_from_file_location("amrtest_" + Path(name).stem,
                                                  PROJECT_ROOT / "scripts" / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def model(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("cpss")
    rng = np.random.default_rng(0)
    n = 150
    ids = [f"g{i:03d}" for i in range(n)]
    signal = rng.random(n) < 0.4
    y = np.where(rng.random(n) < 0.92, signal, ~signal).astype(int)
    rows = [("SIGNAL", signal.astype(int))]
    rows += [(f"U{j}", (rng.random(n) < rng.uniform(0.3, 0.7)).astype(int)) for j in range(80)]
    with open(tmp / "u.rtab", "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(ids) + "\n")
        for seq, bits in rows:
            f.write(seq + "\t" + "\t".join(map(str, bits)) + "\n")
    matrix_store.build_store(tmp / "u.rtab", tmp / "store", min_support=5)
    genomes = pd.DataFrame({"Genome ID": ids, "label": y, "lineage": np.arange(n) // 5})
    matrix_store.build_model_matrix(tmp / "store", genomes, tmp / "model", min_support=5)
    mm = matrix_store.ModelMatrix(tmp / "model")
    store = matrix_store.Store(tmp / "store")
    m = mm.members()
    seq = dict(zip(m["unitig_index"], store.sequences(m["unitig_index"]), strict=True))
    pid = {seq[u]: int(p) for u, p in zip(m["unitig_index"], m["pattern_id"], strict=True)}
    r = _script("04_nested_cv.py")
    assert r.run_folds(mm, tmp / "cv", CV)["evaluable"]
    r.run_final(mm, tmp / "cv", HPO, threads=1)
    return mm, pid, tmp


# ---- CPSS ------------------------------------------------------------------
def test_chi2_prefilter_matches_scipy(model):
    mm, _, _ = model
    table = cpss.chi2_prefilter(mm, 10)
    y = mm.labels.astype(bool)
    for row in table.itertuples():
        x = mm.pattern(row.pattern_id).astype(bool)
        obs = [[(x & y).sum(), (~x & y).sum()], [(x & ~y).sum(), (~x & ~y).sum()]]
        assert row.chi2 == pytest.approx(chi2_contingency(obs, correction=False)[0])
    assert table["chi2"].is_monotonic_decreasing
    assert len(cpss.chi2_prefilter(mm, 10 ** 6)) == mm.n_patterns


def test_half_samples_are_disjoint_stratified_and_seeded():
    y = np.array([1] * 21 + [0] * 40)
    a, b = cpss.half_samples(y, 0, 3)
    assert not set(a) & set(b)
    assert (y[a].sum(), (y[a] == 0).sum()) == (10, 20) == (y[b].sum(), (y[b] == 0).sum())
    assert (cpss.half_samples(y, 0, 3)[0] == a).all() and not (cpss.half_samples(y, 0, 4)[0] == a).all()
    assert cpss.rescaled_colsample(0.001, 2_000_000, 5000) == pytest.approx(0.4)
    assert cpss.rescaled_colsample(0.01, 2_000_000, 5000) == 1.0
    assert cpss.mb_bound(50, 0.6, 5000) == pytest.approx(2.5)


def test_cpss_end_to_end(model, tmp_path):
    mm, pid, tmp = model
    cv = tmp_path / "cv"
    shutil.copytree(tmp / "cv", cv)
    # The tiny search may end with a one-tree final model; CPSS reuses whatever the
    # final model chose, so give it a final record with 20 trees.
    rec = json.loads((cv / "final" / "record.json").read_text())
    rec["n_trees"] = 20
    rec["params"].update(colsample_bytree=0.2, max_depth=3, min_child_weight=1.0, gamma=0.0)
    (cv / "final" / "record.json").write_text(json.dumps(rec))
    s13 = _script("13_cpss.py")
    out = tmp_path / "cpss"
    s13.prefilter(mm, out, CPSS)
    s13.run(mm, cv, out, CPSS, chunk=1, threads=1)
    with pytest.raises(FileNotFoundError, match="1 CPSS chunk"):
        s13.select(mm, cv, out, CPSS, top_gain=5, candidates_file=tmp_path / "c.csv",
                   layers_dir=tmp_path / "layers")
    s13.run(mm, cv, out, CPSS, threads=1)
    s = s13.select(mm, cv, out, CPSS, top_gain=5, candidates_file=tmp_path / "c.csv",
                   layers_dir=tmp_path / "layers")
    assert s["cpss_colsample_bytree"] == pytest.approx(min(1.0, 0.2 * mm.n_patterns / 30))
    table = pd.read_csv(out / "cpss.csv").set_index("pattern_id")
    assert table.loc[pid["SIGNAL"], "pi"] == 1.0 and table.loc[pid["SIGNAL"], "stable"]
    assert table["n_selected"].sum() <= 2 * CPSS["n_pairs"] * CPSS["q"]
    cands = pd.read_csv(tmp_path / "c.csv").set_index("pattern_id")
    assert cands.loc[pid["SIGNAL"], "source"] == "both"
    assert set(cands.index) >= set(table.index[table["stable"]])
    layer = pd.read_csv(tmp_path / "layers" / "cpss.csv").set_index("pattern_id")
    assert set(layer.index) == set(cands.index)
    assert (layer["passes"] == (layer["pi"] >= 0.6)).all()
    assert s["mb_bound"] == pytest.approx(25 / (0.2 * 30)) and s["n_fits"] == 12


def test_a_cpss_chunk_of_another_final_model_is_redone(model, tmp_path):
    """A chunk file written for another final model is neither kept nor read."""
    mm, pid, tmp = model
    cv = tmp_path / "cv"
    shutil.copytree(tmp / "cv", cv)
    rec_file = cv / "final" / "record.json"
    rec = json.loads(rec_file.read_text())
    rec["n_trees"] = 5
    rec_file.write_text(json.dumps(rec))
    s13 = _script("13_cpss.py")
    out = tmp_path / "cpss"
    s13.prefilter(mm, out, CPSS)
    s13.run(mm, cv, out, CPSS, threads=1)
    f = cpss.chunk_file(out, 0)
    key, before = str(np.load(f)["key"]), f.stat().st_mtime_ns
    s13.run(mm, cv, out, CPSS, threads=1)
    assert f.stat().st_mtime_ns == before                    # the same final model: kept
    rec["n_trees"] = 6                                        # a new final model
    rec_file.write_text(json.dumps(rec))
    with pytest.raises(ValueError, match="written for another final model"):
        s13.select(mm, cv, out, CPSS, top_gain=5, candidates_file=tmp_path / "c.csv",
                   layers_dir=tmp_path / "layers")
    s13.run(mm, cv, out, CPSS, threads=1)
    assert str(np.load(f)["key"]) != key


# ---- pyseer ----------------------------------------------------------------
def test_rtab_kinship_and_samples(model, tmp_path):
    mm, _, _ = model
    ids = mm.genomes["Genome ID"].tolist()
    bits = mm.columns([0, 2]).T
    pl.write_rtab(tmp_path / "x.Rtab", ["p0", "p2"], ids, bits)
    back = pd.read_csv(tmp_path / "x.Rtab", sep="\t", index_col=0)
    assert list(back.columns) == ids and (back.to_numpy() == bits).all()
    kin, bg = pl.unitig_samples(mm, 10, 3)
    members = mm.members().sort_values("unitig_index")["pattern_id"].to_numpy()
    assert (kin == members[::10]).all() and len(bg) <= 3
    g = mm.columns(kin).astype(int)
    assert (pl.kinship(mm, kin, block=3) == g @ g.T).all()


def test_pyseer_sample_names_keep_their_digits(tmp_path):
    # pyseer reads the sample column with pandas' type inference (index_col=0)
    ids = ["550.2110", "550.211", "550.21"]
    for names, kept in ((ids, False), (pl.sample_names(ids), True)):
        pd.DataFrame({"samples": names, "resistant": [1, 0, 1]}).to_csv(
            tmp_path / "p.tsv", sep="\t", index=False)
        read = pd.read_csv(tmp_path / "p.tsv", sep="\t", index_col=0).index.astype(str)
        assert (read.tolist() == list(names)) is kept


def test_lambda_and_qq():
    p = np.linspace(0.0005, 0.9995, 1000)
    assert pl.genomic_lambda(p) == pytest.approx(1.0, abs=0.01)
    assert pl.genomic_lambda(chi2.sf(2 * chi2.isf(p, 1), 1)) == pytest.approx(2.0, abs=0.01)
    qq = pl.qq_points([0.01, 0.5, np.nan, 0.1])
    assert len(qq) == 3 and qq["observed"].iloc[0] == pytest.approx(2.0)


def _fake_assoc(path, patterns, pvalues):
    pd.DataFrame({"variant": [f"p{p}" for p in patterns], "af": 0.3, "filter-pvalue": pvalues,
                  "lrt-pvalue": pvalues, "beta": 0.2, "beta-std-err": 0.05, "variant_h2": 0.01,
                  "notes": ""}).to_csv(path, sep="\t", index=False)


def test_pyseer_prep_and_post(model, tmp_path, monkeypatch):
    mm, pid, _ = model
    s14 = _script("14_pyseer.py")
    pd.DataFrame({"pattern_id": [pid["SIGNAL"], 3, 4]}).to_csv(tmp_path / "c.csv", index=False)
    pd.DataFrame({"pattern_id": [pid["SIGNAL"], 5, 6, 7]}).to_csv(tmp_path / "pre.csv",
                                                                   index=False)
    cfg = {"kinship_every": 10, "alpha": 0.05, "background_max": 20}
    out = tmp_path / "pyseer"
    info = s14.prep(mm, out, cfg, prefilter_file=tmp_path / "pre.csv",
                    candidates_file=tmp_path / "c.csv", cpu=2)
    tested = sorted({pid["SIGNAL"], 3, 4, 5, 6, 7})
    assert info["n_tested"] == 6
    assert pd.read_csv(out / "tested.Rtab", sep="\t", index_col=0).index.tolist() == [
        f"p{p}" for p in tested]
    k = pd.read_csv(out / "kinship.tsv", sep="\t", index_col=0)
    assert (k.to_numpy() == k.to_numpy().T).all() and list(k.index) == list(k.columns)
    assert list(k.columns) == pl.sample_names(mm.genomes["Genome ID"].astype(str))
    script = (out / "run_pyseer.sh").read_text()
    assert "--lmm" in script and "--min-af 0 --max-af 1" in script and "--cpu 2" in script
    p = {t: 0.5 for t in tested}
    p[pid["SIGNAL"]] = 1e-6                               # below 0.05 / 6
    p[3] = 0.01                                           # above it
    _fake_assoc(out / "tested_assoc.tsv", tested, [p[t] for t in tested])
    _fake_assoc(out / "background_assoc.tsv", [0, 1, 2], [0.2, 0.5, 0.9])
    s = s14.post(out, cfg, candidates_file=tmp_path / "c.csv", layers_dir=tmp_path / "layers")
    layer = pd.read_csv(tmp_path / "layers" / "pyseer.csv").set_index("pattern_id")
    assert layer["passes"].to_dict() == {pid["SIGNAL"]: True, 3: False, 4: False}
    assert s["bonferroni_threshold"] == pytest.approx(0.05 / 6) and s["n_significant"] == 1
    _fake_assoc(out / "tested_assoc.tsv", tested[:-1], [p[t] for t in tested[:-1]])
    with pytest.raises(SystemExit, match="tested 5 of 6"):
        s14.post(out, cfg, candidates_file=tmp_path / "c.csv", layers_dir=tmp_path / "layers")


@pytest.mark.skipif(shutil.which("pyseer") is None, reason="pyseer is not installed here")
def test_real_pyseer(model, tmp_path):
    mm, pid, _ = model
    s14 = _script("14_pyseer.py")
    pd.DataFrame({"pattern_id": [pid["SIGNAL"], 3]}).to_csv(tmp_path / "c.csv", index=False)
    pd.DataFrame({"pattern_id": list(range(10))}).to_csv(tmp_path / "pre.csv", index=False)
    cfg = {"kinship_every": 2, "alpha": 0.05, "background_max": 20}
    out = tmp_path / "pyseer"
    s14.prep(mm, out, cfg, prefilter_file=tmp_path / "pre.csv",
             candidates_file=tmp_path / "c.csv", cpu=1)
    s14.lmm(out)
    s = s14.post(out, cfg, candidates_file=tmp_path / "c.csv", layers_dir=tmp_path / "layers")
    layer = pd.read_csv(tmp_path / "layers" / "pyseer.csv").set_index("pattern_id")
    assert layer.loc[pid["SIGNAL"], "passes"] and s["pyseer_version"]
    assert json.loads((out / "pyseer_summary.json").read_text())["n_tested"] >= 10

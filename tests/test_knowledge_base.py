#!/usr/bin/env python3
"""Knowledge base (lib.knowledge_base, build_kb.py) built from the real outputs of
the steps on a small synthetic organism (the CARD layer and pyseer's association
table are written by hand)."""
import importlib.util
import json
import shutil
import sqlite3
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
from lib import knowledge_base as kb  # noqa: E402
from lib.card_layer import revcomp  # noqa: E402
from lib.matrix_store import sha256_file  # noqa: E402

pytestmark = pytest.mark.unit

ORG, AB = "ecoli", "ciprofloxacin"
MID = f"{ORG}__{AB}"


def _script(name):
    spec = importlib.util.spec_from_file_location("amrtest_" + Path(name).stem,
                                                  PROJECT_ROOT / "scripts" / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _config(root: Path) -> dict:
    t = {k: str(root / v) for k, v in {
        "panel_dir": "panel", "genome_qc_dir": "{organism}/qc", "lineage_dir": "{organism}/lineage",
        "metadata_file": "{organism}/amr_phenotypes.csv", "matrix_dir": "{organism}/{antibiotic}/mm",
        "cv_dir": "{organism}/{antibiotic}/cv", "rgi_dir": "{organism}/rgi",
        "candidates_file": "{organism}/{antibiotic}/candidates.csv",
        "card_layer_dir": "{organism}/{antibiotic}/card", "layers_dir": "{organism}/{antibiotic}/layers",
        "grades_dir": "{organism}/{antibiotic}/grades", "cpss_dir": "{organism}/{antibiotic}/cpss",
        "pyseer_dir": "{organism}/{antibiotic}/pyseer",
        "label_permutation_dir": "{organism}/{antibiotic}/lp", "cross_model_dir": "cross",
        "external_dir": "{organism}/external", "context_dir": "{organism}/context",
        "unitig_store_dir": "{organism}/store", "databases_manifest": "databases.json",
        "kb_dir": "kb", "reports_dir": "reports",
        "model_reports_dir": "reports/{organism}/{antibiotic}"}.items()}
    text = PROJECT_ROOT / "docs" / "V1_PROTOKOL.md"
    return {
        "paths_organism": t,
        "protocol": {"version": "1.0", "sha256": sha256_file(text) if text.exists() else "x"},
        "panel": {"min_minority": 20}, "unitig": {"k": 31, "min_support": 5},
        "hpo": {"n_trials": 3, "search_max_genomes": 80, "pruner_startup_trials": 1,
                "pruner_warmup_rounds": 5, "early_stopping_rounds": 10, "max_rounds": 100,
                "batch_rows": 50},
        "cv": {"n_repeats": 1, "n_folds": 3, "min_minority_per_test_fold": 5,
               "max_seed_attempts": 20, "threshold": 0.5, "n_bootstrap": 30,
               "reliability_bins": 5},
        "card": {"near_universal": 0.95, "max_located_genomes": 3, "min_overlap": 0.5},
        "prevalence": {"min_delta": 0.10, "alpha": 0.05},
        "mda": {"n_permutations": 20, "seed": 0, "alpha": 0.05, "cluster_r": 0.9},
        "cpss": {"prefilter": 20, "n_pairs": 4, "q": 5, "pi_threshold": 0.6, "seed": 0,
                 "chunk": 2},
        "candidates": {"top_gain": 5},
        "pyseer": {"kinship_every": 5, "alpha": 0.05, "background_max": 10},
        "label_permutation": {"n_permutations": 6, "seed": 0, "chunk": 3, "alpha": 0.05},
        "grading": {"rule": "allele_aware"},
    }


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    root = tmp_path_factory.mktemp("kb")
    config = _config(root)
    p = {k: Path(v.format(organism=ORG, antibiotic=AB))
         for k, v in config["paths_organism"].items()}
    rng = np.random.default_rng(0)
    n = 150
    ids = [f"562.{i:04d}" for i in range(n)]
    signal = rng.random(n) < 0.4
    y = np.where(rng.random(n) < 0.92, signal, ~signal).astype(int)
    bases = np.array(list("ACGT"))
    seqs = {"SIGNAL": "".join(rng.choice(bases, 40))}
    rows = [(seqs["SIGNAL"], signal.astype(int))]
    for j in range(40):
        seqs[f"U{j}"] = "".join(rng.choice(bases, 35))
        rows.append((seqs[f"U{j}"], (rng.random(n) < rng.uniform(0.3, 0.7)).astype(int)))
    # inputs of the panel: phenotypes, QC, lineages, decisions
    for d in ("genome_qc_dir", "lineage_dir", "panel_dir", "rgi_dir"):
        p[d].mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"Genome ID": ids, AB: y, "ampicillin": np.where(np.arange(n) < 20, 1, np.nan)}
                 ).to_csv(p["metadata_file"], index=False)
    # the data snapshot (00a): NCBI identifiers and snapshot.json
    meta = p["metadata_file"].parent
    pd.DataFrame({"genome_id": ids, "genome_name": "Escherichia coli",
                  "taxon_id": ["562"] * (n - 1) + ["83334"], "genome_status": "WGS",
                  "assembly_accession": ["GCA_1.1"] + [""] * (n - 1),
                  "sra_accession": ["SRR1,SRR2"] + [""] * (n - 1),
                  "biosample_accession": [f"SAMN{i}" for i in range(n)],
                  "bioproject_accession": "", "contigs": 80, "genome_length": 5_000_000}).to_csv(
        meta / "genomes.csv", index=False)
    (meta / "snapshot.json").write_text(json.dumps({
        "frozen_at": "2026-10-20T10:00:00+00:00", "n_genomes": n,
        "query": {"source": "https://www.bv-brc.org/api", "api_version": "1.9.3",
                  "queried_at": "2026-10-20T08:00:00+00:00", "genome_query": "g",
                  "record_query": "r"}}))
    pd.DataFrame({"genome_id": ids, "completeness": 99.0, "contamination": 0.5, "n50": 90000,
                  "n_contigs": 80, "total_length": 5_000_000, "pass_completeness": True,
                  "pass_contamination": True, "pass_n50": True, "pass_contigs": True,
                  "pass_overall": True}).to_csv(p["genome_qc_dir"] / f"02d_genome_qc_{ORG}.csv",
                                                index=False)
    pd.DataFrame({"Genome ID": ids, "Cluster": np.arange(n) // 5}).to_csv(
        p["lineage_dir"] / "poppunk_clusters.csv", index=False)
    pd.DataFrame([{"organism": ORG, "antibiotic": AB, "drug_class": "quinolones", "n_tested": n,
                   "n_eligible": n, "n_resistant": int(y.sum()), "n_susceptible": int(n - y.sum()),
                   "minority": int(min(y.sum(), n - y.sum())), "decision": "included", "reason": ""},
                  {"organism": ORG, "antibiotic": "ampicillin", "drug_class": "penicillins",
                   "n_tested": 20, "n_eligible": 20, "n_resistant": 20, "n_susceptible": 0,
                   "minority": 0, "decision": "excluded", "reason": "minority below 20"}]
                 ).to_csv(p["panel_dir"] / "panel_decisions.csv", index=False)
    # 03u, 04
    with open(root / "u.rtab", "w") as f:
        f.write("Unitig_sequence\t" + "\t".join(ids) + "\n")
        for seq, bits in rows:
            f.write(seq + "\t" + "\t".join(map(str, bits)) + "\n")
    (root / "versions.json").write_text(json.dumps({"unitig_caller": "1.3.1"}))
    matrix_store.build_store(root / "u.rtab", p["unitig_store_dir"], min_support=5)
    genomes = pd.DataFrame({"Genome ID": ids, "label": y, "lineage": np.arange(n) // 5})
    matrix_store.build_model_matrix(p["unitig_store_dir"], genomes, p["matrix_dir"],
                                    min_support=5)
    mm = matrix_store.ModelMatrix(p["matrix_dir"])
    r04 = _script("04_nested_cv.py")
    assert r04.run_folds(mm, p["cv_dir"], config["cv"])["evaluable"]
    for arm in folds.ARMS:
        for k in range(3):
            r04.run_unit(mm, p["cv_dir"], arm, 1, k, config["hpo"], threads=1)
    r04.run_final(mm, p["cv_dir"], config["hpo"], threads=1)
    r04.run_metrics(mm, p["cv_dir"], config["cv"])
    # 13 (candidates), 10, 12, 12b
    s13 = _script("13_cpss.py")
    s13.prefilter(mm, p["cpss_dir"], config["cpss"])
    s13.run(mm, p["cv_dir"], p["cpss_dir"], config["cpss"], threads=1)
    s13.select(mm, p["cv_dir"], p["cpss_dir"], config["cpss"], top_gain=5,
               candidates_file=p["candidates_file"], layers_dir=p["layers_dir"])
    argv, saved = ["x", "--organism", ORG, "--antibiotic", AB], sys.argv
    sys.argv = argv
    try:
        for name in ("10_prevalence.py", "12_mda.py"):
            m = _script(name)
            m.load_config = lambda: config
            m.main()
    finally:
        sys.argv = saved
    s12b = _script("12b_label_permutation.py")
    s12b.run(mm, p["cv_dir"], p["label_permutation_dir"], config["label_permutation"], threads=1)
    s12b.metrics(mm, p["cv_dir"], p["label_permutation_dir"], config["label_permutation"])
    s12b.across(config)
    # 14 with a hand-written association table
    s14 = _script("14_pyseer.py")
    s14.prep(mm, p["pyseer_dir"], config["pyseer"], prefilter_file=p["cpss_dir"] / "prefilter.csv",
             candidates_file=p["candidates_file"], cpu=1)
    tested = pd.read_csv(p["pyseer_dir"] / "tested_patterns.csv")["pattern_id"]
    for name, pats in (("tested", tested), ("background", tested[:5])):
        pd.DataFrame({"variant": [f"p{x}" for x in pats], "af": 0.3, "filter-pvalue": 0.5,
                      "lrt-pvalue": np.linspace(1e-6, 0.9, len(pats)), "beta": 0.1,
                      "beta-std-err": 0.05, "variant_h2": 0.01, "notes": ""}).to_csv(
            p["pyseer_dir"] / f"{name}_assoc.tsv", sep="\t", index=False)
    s14.post(p["pyseer_dir"], config["pyseer"], candidates_file=p["candidates_file"],
             layers_dir=p["layers_dir"])
    # 08/09 by hand: the SIGNAL unitig lies in gyrA (a variant hit carrying the allele)
    cands = pd.read_csv(p["candidates_file"])["pattern_id"].astype(int)
    store = matrix_store.Store(p["unitig_store_dir"])
    members = mm.members()
    members = members[members["pattern_id"].isin(set(cands))]
    seq = store.sequences(members["unitig_index"])
    is_signal = [s == seqs["SIGNAL"] for s in seq]
    st = np.where(is_signal, "b", "no_card_hit")
    card = pd.DataFrame({"pattern_id": members["pattern_id"], "unitig_index": members["unitig_index"],
                         "sequence": seq, "length": [len(s) for s in seq],
                         "located_genomes": ids[0]})
    for mode, state in (("allele_aware", st), ("homolog_only",
                                                np.where(is_signal, "card_hit_without_b",
                                                         "no_card_hit"))):
        card[f"{mode}_state"] = state
        card[f"{mode}_reasons"] = np.where(state == "card_hit_without_b", "variant_not_counted", "")
        card[f"{mode}_aros"] = np.where(is_signal, "3003294", "")
        card[f"{mode}_n_b"] = (state == "b").astype(int) * 3
    p["card_layer_dir"].mkdir(parents=True)
    card.to_csv(p["card_layer_dir"] / "card_unitigs.csv", index=False)
    pats = card.groupby("pattern_id").agg(n_members=("unitig_index", "size")).reset_index()
    for mode in ("allele_aware", "homolog_only"):
        best = card.groupby("pattern_id")[f"{mode}_state"].agg(
            lambda s: "b" if (s == "b").any() else ("card_hit_without_b"
                                                    if (s == "card_hit_without_b").any()
                                                    else "no_card_hit"))
        pats[f"{mode}_state"] = pats["pattern_id"].map(best)
        pats[f"{mode}_reasons"] = np.where(pats[f"{mode}_state"] == "card_hit_without_b",
                                           "variant_not_counted", "")
        pats[f"{mode}_aros"] = ""
    pats.to_csv(p["card_layer_dir"] / "card_patterns.csv", index=False)
    pd.DataFrame({"genome_id": ids[0], "aro": ["3003294"], "aro_name": ["Escherichia coli gyrA"],
                  "model_type": ["variant"], "gene_family": ["fluoroquinolone resistant gyrA"],
                  "drug_class": ["fluoroquinolone antibiotic"],
                  "mechanism": ["antibiotic target alteration"]}).to_csv(
        p["rgi_dir"] / "rgi_hits.csv", index=False)
    (p["rgi_dir"] / "rgi_summary.json").write_text(json.dumps(
        {"rgi_version": "6.0.8", "card_version": "4.0.1"}))
    m14b = _script("14b_grading.py")
    m14b.load_config = lambda: config
    sys.argv = argv
    try:
        m14b.main()
    finally:
        sys.argv = saved
    # 16 by hand: one AMRFinderPlus call and the model's row of the comparison
    from lib import external as ex
    ext = p["external_dir"]
    ext.mkdir(parents=True)
    pd.DataFrame([(ids[0], "gyrA_S83L", "AMR", "POINT", "core", "QUINOLONE", "QUINOLONE")],
                 columns=ex.CALL_COLUMNS).to_csv(ext / "amrfinder_calls.csv", index=False)
    pd.DataFrame({"genome_id": ids}).to_csv(ext / "amrfinder_genomes.csv", index=False)
    oof = pd.read_csv(p["cv_dir"] / "oof_predictions.csv", dtype={"genome_id": str})
    ex.compare_model(MID, pd.Series(y, index=ids), {"model": ex.predict_model(oof, ids, 0.5)}
                     ).to_csv(ext / "comparison.csv", index=False)
    (ext / "versions.json").write_text(json.dumps({"amrfinderplus": "4.2.7",
                                                   "amrfinderplus_database": "2025-07-16.1",
                                                   "resfinder": "4.5.0"}))
    # 18 by hand: the SIGNAL unitig was searched and sits on a plasmid in half of its hits
    p["context_dir"].mkdir(parents=True)
    pd.DataFrame([{"unitig_id": kb.unitig_id(seqs["SIGNAL"])[0], "organism_id": ORG,
                   "source": "ncbi_nt_remote", "queried_on": "2026-10-03", "nt_release": "x",
                   "n_hits": 2, "best_accession": "CP000001", "best_title": "E. coli plasmid",
                   "best_identity": 100.0, "best_coverage": 1.0, "best_evalue": 1e-20,
                   "gene": "gyrA", "product": "DNA gyrase subunit A", "plasmid_share": 0.5}]
                 ).to_csv(p["context_dir"] / "unitig_context.csv", index=False)
    (p["lineage_dir"] / "versions.json").write_text(json.dumps({"poppunk": "2.7.8",
                                                                 "graph_tool": "2.98"}))
    (p["genome_qc_dir"] / "versions.json").write_text(json.dumps({"checkm2": "1.1.0",
                                                                   "quast": "QUAST v5.3.0"}))
    (root / "card").mkdir()
    (root / "card" / "card.json").write_text(json.dumps({"_version": "4.0.1"}))
    from lib import databases
    databases.record(p["databases_manifest"], "card", root / "card", downloaded_on="2026-11-02")
    out = p["kb_dir"] / "kanit.sqlite"
    report = _script("build_kb.py").build(config, out, kb_version="1.0.0-test")
    return config, p, mm, seqs, report, out


def test_unitig_id_is_stable_and_orientation_free():
    i, c = kb.unitig_id("acgtTT")
    assert kb.unitig_id(revcomp("ACGTTT")) == (i, c) and c == "AAACGT"
    assert i.startswith("KU") and len(i) == 22


def test_build_loads_every_layer_and_rechecks_the_grades(built):
    _, p, mm, _, report, out = built
    n_cand = len(pd.read_csv(p["candidates_file"]))
    t = report["tables"]
    assert t["model"] == 1 and t["candidate"] == n_cand and t["grade"] == 2 * n_cand
    assert report["grades_rechecked"] == 2 * n_cand and report["sha256"] == sha256_file(out)
    assert t["panel_decision"] == 2 and t["genome"] == 150 and t["pyseer_result"] >= n_cand
    with sqlite3.connect(out) as conn:
        assert conn.execute("PRAGMA foreign_key_check").fetchall() == []
        rel = conn.execute("SELECT protocol_version, card_version, tools FROM release").fetchone()
        assert rel[0] == "1.0" and rel[1] == "4.0.1" and "pyseer" in json.loads(rel[2])
        assert conn.execute("SELECT count(*) FROM v_biomarker").fetchone()[0] == n_cand
        auc = conn.execute("SELECT roc_auc_lineage_aware, flag FROM v_model").fetchone()
        assert auc[0] > 0.7 and auc[1] in ("", "permutation_not_significant")
        assert conn.execute("SELECT count(*) FROM source_file").fetchone()[0] > 15
        ba = conn.execute("SELECT assessable, balanced_accuracy FROM external_comparison "
                          "WHERE tool = 'model'").fetchone()
        assert ba[0] == 1 and ba[1] > 0.7
        assert conn.execute("SELECT count(*) FROM external_comparison WHERE assessable = 0"
                            ).fetchone()[0] == 4
        assert json.loads(rel[2])["resfinder"] == "4.5.0"
        assert conn.execute("SELECT gene FROM unitig_context").fetchone()[0] == "gyrA"
        tools = json.loads(rel[2])
        assert tools["graph_tool"] == "2.98" and tools["unitig_caller"] == "1.3.1"
        assert tools["rgi"] == "6.0.8" and "xgboost" in tools
        assert (tools["checkm2"], tools["quast"]) == ("1.1.0", "QUAST v5.3.0")
        assert conn.execute("SELECT version, downloaded_on FROM reference_database").fetchone() \
            == ("4.0.1", "2026-11-02")
        assert conn.execute("SELECT ncbi_taxid, assembly_accession, sra_accession, "
                            "biosample_accession FROM genome ORDER BY genome_id").fetchone() \
            == (562, "GCA_1.1", "SRR1,SRR2", "SAMN0")
        assert conn.execute("SELECT count(*) FROM genome WHERE assembly_accession IS NULL "
                            "AND ncbi_taxid = 83334").fetchone()[0] == 1
        snap = conn.execute("SELECT api_version, n_genomes, sha256 FROM data_snapshot").fetchone()
        assert snap[:2] == ("1.9.3", 150) and snap[2] == sha256_file(
            p["metadata_file"].parent / "snapshot.json")
        assert conn.execute("SELECT ncbi_taxids FROM organism").fetchone()[0] == "562"


def test_signal_is_graded_from_its_card_variant_hit(built):
    _, p, mm, seqs, _, out = built
    uid, _ = kb.unitig_id(seqs["SIGNAL"])
    with sqlite3.connect(out) as conn:
        pid, grade_aa, grade_ho, genes = conn.execute(
            "SELECT b.pattern_id, b.grade_allele_aware, b.grade_homolog_only, b.card_genes "
            "FROM v_biomarker b JOIN pattern_member pm USING (model_id, pattern_id) "
            "WHERE pm.unitig_id = ?", (uid,)).fetchone()
        assert genes == "Escherichia coli gyrA"
        assert grade_aa in ("confirmed", "candidate") and grade_ho != "confirmed"
        blob = conn.execute("SELECT carriers FROM pattern WHERE model_id = ? AND pattern_id = ?",
                            (MID, pid)).fetchone()[0]
    bits = np.unpackbits(np.frombuffer(blob, dtype=np.uint8))[:mm.n_genomes]
    assert (bits == mm.pattern(pid)).all()


def test_a_grade_that_does_not_follow_its_layers_stops_the_build(built, tmp_path):
    config, p, _, _, _, out = built
    before = sha256_file(out)
    g = p["grades_dir"] / "grades_patterns.csv"
    original = g.read_text()
    d = pd.read_csv(g, keep_default_na=False)
    d.loc[0, "allele_aware_grade"] = "strong_novel" if d.loc[0, "allele_aware_grade"] != \
        "strong_novel" else "none"
    d.to_csv(g, index=False)
    try:
        with pytest.raises(ValueError, match="does not follow its layers"):
            _script("build_kb.py").build(config, out, kb_version="bad")
    finally:
        g.write_text(original)
    assert sha256_file(out) == before                      # the old file is kept
    shutil.move(p["layers_dir"] / "mda.csv", tmp_path / "mda.csv")
    try:
        with pytest.raises(FileNotFoundError, match="mda.csv"):
            _script("build_kb.py").build(config, out, kb_version="bad")
    finally:
        shutil.move(tmp_path / "mda.csv", p["layers_dir"] / "mda.csv")


def test_outputs_follow_the_contract(built):
    """Every table the steps wrote here passes its schema in config/output_contract.yaml."""
    from lib import contract
    config = built[0]
    checked = []
    for tid, t in contract.load()["tables"].items():
        try:
            path = contract.table_path(tid, config, ORG, AB)
        except KeyError:                      # a location this synthetic run does not set
            continue
        if path.exists():
            assert contract.validate_csv(path, t) == [], (tid, path)
            checked.append(tid)
    assert {"panel_decisions", "genome_qc", "model_patterns", "oof_predictions", "candidates",
            "cpss_layer", "prevalence_layer", "mda_layer", "pyseer_layer", "card_unitigs",
            "grades_patterns", "rgi_hits", "external_comparison"} <= set(checked), checked


def test_reports_from_the_tables(built, monkeypatch):
    """Figures, the model report and the numbers file come from the steps' CSVs only."""
    from lib import contract
    config = built[0]
    m = _script("reports.py")
    monkeypatch.setattr(m, "load_config", lambda: config)
    argv, saved = ["x", "--entry", "main", "--organisms", ORG], sys.argv
    sys.argv = argv
    try:
        m.main()
    finally:
        sys.argv = saved
    rep = Path(config["paths_organism"]["reports_dir"])
    model_dir = rep / ORG / AB
    names = {f.stem for f in model_dir.glob("*.png")}
    assert {"matrix", "cross_validation", "cpss", "prevalence", "mda", "card", "grading"} <= names
    page = (model_dir / "report.html").read_text()
    for section in ("Cross-validation", "Evidence layers", "Best candidates", "base64,"):
        assert section in page
    assert (rep / "panel.png").exists() and (rep / ORG / "genome_qc.png").exists()
    nums = rep / "tez_sayilari.csv"
    assert contract.validate_csv(nums, contract.load()["tables"]["thesis_numbers"]) == []
    t = pd.read_csv(nums).set_index("key")
    assert int(t.loc[f"{MID}.n_genomes", "value"]) == 150
    assert t.loc[f"{MID}.roc_auc.lineage_aware", "source_table"] == "repeat_metrics"

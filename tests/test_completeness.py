#!/usr/bin/env python3
"""Completeness of a workflow run (scripts/completeness.py)."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
pytestmark = pytest.mark.unit


def _script():
    spec = importlib.util.spec_from_file_location("amrtest_completeness",
                                                  PROJECT_ROOT / "scripts" / "completeness.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def _write(m, c, p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)
    if c.kind == "json":
        p.write_text("{}")
    elif c.kind == "csv":
        cols = list(c.columns) or ["x"]
        pd.DataFrame([dict.fromkeys(cols, 1)] * max(c.min_rows, 1)).to_csv(p, index=False)
    else:
        p.write_text("x")


def test_run_checks(tmp_path):
    m = _script()
    keys = {"metadata_file": "{organism}/meta/amr_phenotypes.csv", "genome_qc_dir": "{organism}/qc",
            "lineage_dir": "{organism}/lin", "unitig_store_dir": "{organism}/store",
            "rgi_dir": "{organism}/rgi", "external_dir": "{organism}/ext",
            "matrix_dir": "{organism}/{antibiotic}/mm", "cv_dir": "{organism}/{antibiotic}/cv",
            "label_permutation_dir": "{organism}/{antibiotic}/lp",
            "cpss_dir": "{organism}/{antibiotic}/cpss",
            "candidates_file": "{organism}/{antibiotic}/candidates.csv",
            "layers_dir": "{organism}/{antibiotic}/layers",
            "pyseer_dir": "{organism}/{antibiotic}/pyseer",
            "card_layer_dir": "{organism}/{antibiotic}/card",
            "grades_dir": "{organism}/{antibiotic}/grades", "context_dir": "{organism}/ctx",
            "panel_dir": "panel", "cross_model_dir": "cross", "kb_dir": "kb"}
    config = {"paths_organism": {k: str(tmp_path / v) for k, v in keys.items()}}
    (tmp_path / "panel").mkdir()
    pd.DataFrame({"organism": ["ecoli", "ecoli", "ecoli"],
                  "antibiotic": ["ampicillin", "colistin", "tetracycline"],
                  "decision": ["included", "included", "excluded"]}).to_csv(
        tmp_path / "panel" / "panel_decisions.csv", index=False)
    for c in m.ORGANISM + m.MODELLED_ORGANISM:
        _write(m, c, m.path_of(c, config, "ecoli"))
    for ab, evaluable in (("ampicillin", True), ("colistin", False)):
        for c in m.MODEL:
            if evaluable or not c.evaluable_only:
                _write(m, c, m.path_of(c, config, "ecoli", ab))
        m.path_of(m.MODEL[3], config, "ecoli", ab).write_text(json.dumps({"evaluable": evaluable}))
    for c in m.ACROSS:
        _write(m, c, m.path_of(c, config))
    t = m.run_checks("main", ["ecoli"], config)
    assert set(t["status"]) == {"ok", "n/a"} and len(t[t["unit"] == "ecoli__tetracycline"]) == 0
    assert (t.loc[t["unit"] == "ecoli__colistin", "status"] == "n/a").sum() == sum(
        c.evaluable_only for c in m.MODEL)
    # a missing and an invalid output
    (tmp_path / "ecoli" / "ampicillin" / "grades" / "grades_summary.json").unlink()
    pd.DataFrame({"genome_id": [1]}).to_csv(tmp_path / "ecoli" / "ampicillin" / "cv" /
                                            "oof_predictions.csv", index=False)
    t = m.run_checks("main", ["ecoli"], config).set_index(["unit", "step", "status"])
    bad = t[t.index.get_level_values("status").isin(["missing", "invalid"])]
    assert sorted(bad.index.droplevel("unit")) == [("cv", "invalid"), ("grading", "missing")]
    assert "no column ['p']" in bad["detail"].tolist()
    mat = m.matrix(m.run_checks("main", ["ecoli"], config)).set_index("unit")
    assert mat.loc["ecoli__ampicillin", "grading"] == "missing"
    assert mat.loc["ecoli__colistin", "grading"] == "n/a"
    # the snapshot alone for the DOWNLOAD entry, and the gate
    assert set(m.run_checks("DOWNLOAD", ["ecoli"], config)["step"]) == {"snapshot"}
    out = tmp_path / "run"
    out.mkdir()
    (out / "completeness_summary.json").write_text(json.dumps({
        "complete": False, "percent_ok": 50.0, "applicable": 2,
        "failures": [{"unit": "u", "step": "s", "file": "f", "status": "missing", "detail": ""}]}))
    r = subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "completeness.py"), "gate",
                        "--out-dir", str(out)], capture_output=True, text=True)
    assert r.returncode == 1 and "INCOMPLETE RUN" in r.stderr

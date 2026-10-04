#!/usr/bin/env python3
"""Completeness of a workflow run (scripts/completeness.py) and the output contract
(config/output_contract.yaml, lib/contract.py)."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import contract  # noqa: E402
from lib.config import load_config  # noqa: E402

pytestmark = pytest.mark.unit


def _script():
    spec = importlib.util.spec_from_file_location("amrtest_completeness",
                                                  PROJECT_ROOT / "scripts" / "completeness.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def _value(f):
    if "enum" in f:
        return f["enum"][0]
    return {"integer": f.get("minimum", 1), "number": f.get("minimum", 0.5),
            "boolean": True}.get(f.get("type", "string"), "x")


def _write(kind, path: Path, table=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if kind == "csv":
        pd.DataFrame([{f["name"]: _value(f) for f in table["fields"]}]).to_csv(path, index=False)
    elif kind == "json":
        path.write_text("{}")
    else:
        path.write_text("x")


def test_contract_is_consistent():
    c = contract.load()
    used = set()
    for name, step in c["steps"].items():
        assert step["entry"] in ("main", "DOWNLOAD", "CONTEXT", "KB"), name
        assert step["level"] in ("organism", "modelled_organism", "model", "global"), name
        tables = step.get("tables") or []
        assert step["main"] is None or step["main"] in tables, name
        used |= set(tables)
    assert used == set(c["tables"])                  # every schema belongs to a step
    for tid, t in c["tables"].items():
        names = [f["name"] for f in t["fields"]]
        assert len(names) == len(set(names)) and all(f.get("description") for f in t["fields"])
        assert set(t.get("primary_key", [])) <= set(names), tid
        assert all(f.get("type", "string") in ("string", "integer", "number", "boolean")
                   for f in t["fields"]), tid


def test_validation_rules():
    t = {"fields": [{"name": "id", "type": "integer", "required": True},
                    {"name": "p", "type": "number", "minimum": 0, "maximum": 1},
                    {"name": "y", "type": "integer", "enum": [0, 1]},
                    {"name": "ok", "type": "boolean"}], "primary_key": ["id"], "min_rows": 1}
    good = pd.DataFrame({"id": ["1", "2.0"], "p": ["0.5", ""], "y": ["1", "0"],
                         "ok": ["True", "false"]})
    assert contract.validate_frame(good, t) == []
    bad = pd.DataFrame({"id": ["1", "1"], "p": ["1.5", "x"], "y": ["2", "0"], "ok": ["yes", ""],
                        "z": ["?", "?"]})
    problems = " | ".join(contract.validate_frame(bad, t))
    for text in ("undocumented column ['z']", "p: '1.5' above 1", "y: '2' not in",
                 "ok: 'yes' not a boolean", "primary key ['id'] not unique"):
        assert text in problems, text


def test_run_checks(tmp_path):
    m = _script()
    c = contract.load()
    config = {"paths_organism": {k: str(tmp_path / v)
                                 for k, v in load_config()["paths_organism"].items()}}
    panel_path = contract.table_path("panel_decisions", config)
    panel_path.parent.mkdir(parents=True)
    base = {f["name"]: _value(f) for f in c["tables"]["panel_decisions"]["fields"]}
    pd.DataFrame([{**base, "organism": "ecoli", "antibiotic": ab, "decision": d} for ab, d in (
        ("ampicillin", "included"), ("colistin", "included"), ("tetracycline", "excluded"))]
                 ).to_csv(panel_path, index=False)
    for step in c["steps"].values():
        if step["entry"] != "main":
            continue
        units = {"global": [(None, None)], "organism": [("ecoli", None)],
                 "modelled_organism": [("ecoli", None)],
                 "model": [("ecoli", "ampicillin"), ("ecoli", "colistin")]}[step["level"]]
        for org, ab in units:
            for x in m.step_checks(step):
                if x["evaluable_only"] and ab == "colistin":
                    continue
                path = contract.output_path(x["location"], x["file"], config, org, ab)
                if path == panel_path:
                    continue
                _write(x["kind"], path, c["tables"][x["table"]] if x["table"] else None)
    for ab, ev in (("ampicillin", True), ("colistin", False)):
        design = contract.output_path("cv_dir", "cv_design.json", config, "ecoli", ab)
        design.write_text(json.dumps({"evaluable": ev}))
    t = m.run_checks("main", ["ecoli"], config)
    assert set(t["status"]) == {"ok", "n/a"}, t[~t["status"].isin(["ok", "n/a"])].to_dict("records")
    assert "ecoli__tetracycline" not in set(t["unit"])
    assert set(t.loc[t["unit"] == "ecoli__colistin", "step"]) >= {"model_matrix", "cv", "grading"}
    # a missing output and one that breaks its schema
    contract.output_path("grades_dir", "grades_summary.json", config, "ecoli", "ampicillin").unlink()
    oof = contract.table_path("oof_predictions", config, "ecoli", "ampicillin")
    pd.read_csv(oof).assign(p=7).to_csv(oof, index=False)
    t = m.run_checks("main", ["ecoli"], config)
    bad = t[t["status"].isin(["missing", "invalid"])]
    assert sorted(zip(bad["step"], bad["status"], strict=True)) == [("cv", "invalid"),
                                                                    ("grading", "missing")]
    assert "above 1" in bad["detail"].str.cat()
    mat = m.matrix(t).set_index("unit")
    assert mat.loc["ecoli__ampicillin", "grading"] == "missing"
    assert mat.loc["ecoli__colistin", "grading"] == "n/a"
    assert set(m.run_checks("DOWNLOAD", ["ecoli"], config)["step"]) == {"snapshot"}
    out = tmp_path / "run"
    out.mkdir()
    (out / "completeness_summary.json").write_text(json.dumps({
        "complete": False, "percent_ok": 50.0, "applicable": 2,
        "failures": [{"unit": "u", "step": "s", "file": "f", "status": "missing", "detail": ""}]}))
    r = subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "completeness.py"), "gate",
                        "--out-dir", str(out)], capture_output=True, text=True)
    assert r.returncode == 1 and "INCOMPLETE RUN" in r.stderr

#!/usr/bin/env python3
"""Reference database manifest (lib.databases) and versions read from the tools."""
import json
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import databases as db  # noqa: E402
from lib import run_metadata  # noqa: E402

pytestmark = pytest.mark.unit


def _dbs(root: Path) -> dict:
    (root / "card").mkdir()
    (root / "card" / "card.json").write_text(json.dumps({"_version": "4.0.1", "1": {}}))
    (root / "amrfinder").mkdir()
    (root / "amrfinder" / "version.txt").write_text("2025-07-16.1\n")
    (root / "amrfinder" / "AMR.LIB").write_text("x")
    (root / "resfinder").mkdir()
    (root / "resfinder" / "VERSION").write_text("2.4.0\n")
    (root / "checkm2").mkdir()
    (root / "checkm2" / "uniref100.KO.1.dmnd").write_text("y")
    return {"card": root / "card", "amrfinderplus": root / "amrfinder",
            "resfinder": root / "resfinder", "checkm2": root / "checkm2"}


def test_versions_are_read_from_the_databases(tmp_path):
    p = _dbs(tmp_path)
    assert db.detect_version("card", p["card"]) == "4.0.1"
    assert db.detect_version("amrfinderplus", p["amrfinderplus"]) == "2025-07-16.1"
    assert db.detect_version("resfinder", p["resfinder"]) == "2.4.0"
    assert db.detect_version("checkm2", p["checkm2"]) == "uniref100.KO.1.dmnd"


def test_record_verify_and_same_day(tmp_path):
    p = _dbs(tmp_path)
    m = tmp_path / "databases.json"
    for name, path in p.items():
        db.record(m, name, path, source="https://example.org", downloaded_on="2026-11-02")
    doc = db.load(m)
    assert doc["card"]["version"] == "4.0.1" and len(doc["amrfinderplus"]["files"]) == 2
    assert db.verify(m) == [] and db.same_day(m)
    (p["amrfinderplus"] / "AMR.LIB").write_text("changed")
    (p["resfinder"] / "new.fsa").write_text("z")
    problems = db.verify(m)
    assert any(x.startswith("amrfinderplus") for x in problems)
    assert any(x.startswith("resfinder") for x in problems)
    db.record(m, "pointfinder", p["resfinder"], version="v1", downloaded_on="2026-11-03")
    assert not db.same_day(m)
    with pytest.raises(ValueError, match="unknown database"):
        db.record(m, "blast", p["card"])
    (tmp_path / "empty").mkdir()
    with pytest.raises(ValueError, match="no version"):
        db.record(m, "amrfinderplus", tmp_path / "empty")


def test_card_version_comes_from_card_json_not_config(tmp_path):
    (tmp_path / "card.json").write_text(json.dumps({"_version": "4.0.2"}))
    v = run_metadata.collect_versions(config={"card": {"card_json": str(tmp_path / "card.json")},
                                              "provenance": {"card_version": "9.9"}})
    assert v["card_version"] == "4.0.2"
    import yaml
    cfg = yaml.safe_load((PROJECT_ROOT / "config" / "config.yaml").read_text())
    assert "card_version" not in (cfg.get("provenance") or {})

#!/usr/bin/env python3
"""Release package of the knowledge base (kb_package.py) from a small knowledge base
written by hand, and the data dictionary against the schema."""
import base64
import copy
import csv
import hashlib
import importlib.util
import io
import json
import sqlite3
import sys
import zipfile
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import knowledge_base as kb  # noqa: E402
from lib.matrix_store import sha256_file  # noqa: E402

pytestmark = pytest.mark.unit

MID = "ecoli__ciprofloxacin"
CARRIERS = bytes([0b10100000])
UID, SEQ = kb.unitig_id("ACGTTGCAAGT")


def _module():
    spec = importlib.util.spec_from_file_location("amrtest_kb_package",
                                                  PROJECT_ROOT / "scripts" / "kb_package.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


m = _module()


def _dictionary() -> dict:
    return yaml.safe_load(m.DICTIONARY_FILE.read_text(encoding="utf-8"))


def _kb(path: Path, license=None, doi=None) -> Path:
    conn = sqlite3.connect(path)
    kb.create(conn)
    with conn:
        conn.execute("INSERT INTO organism VALUES ('ecoli', 'Escherichia coli', '562', "
                     "'negative', 'Pseudomonadota')")
        conn.execute("INSERT INTO antibiotic VALUES ('ciprofloxacin', 'fluoroquinolone', "
                     "'fluoroquinolone antibiotic', 'Watch')")
        conn.executemany("INSERT INTO genome (genome_id, organism_id, qc_pass, ncbi_taxid, "
                         "lineage_cluster) VALUES (?, 'ecoli', 1, 562, ?)",
                         [("562.1", "1"), ("562.2", "1"), ("562.3", None)])
        conn.execute("INSERT INTO model VALUES (?, 'ecoli', 'ciprofloxacin', 3, 2, 1, 2, 0.667, "
                     "40, 12, 1, 150, '{\"max_depth\": 4}')", (MID,))
        conn.executemany("INSERT INTO model_genome VALUES (?, ?, ?, ?, ?)",
                         [(MID, "562.1", 0, 1, "1"), (MID, "562.2", 1, 0, "1"),
                          (MID, "562.3", 2, 1, "2")])
        conn.execute("INSERT INTO label_permutation VALUES (?, 0.91, 0.5, NULL, NULL, 100, "
                     "0.0099, 0.02, '')", (MID,))
        conn.execute("INSERT INTO unitig VALUES (?, ?, ?)", (UID, SEQ, len(SEQ)))
        conn.execute("INSERT INTO pattern VALUES (?, 7, 1, 2, ?)", (MID, CARRIERS))
        conn.execute("INSERT INTO pattern VALUES (?, 8, 3, 1, NULL)", (MID,))
        conn.execute("INSERT INTO pattern_member VALUES (?, 7, ?)", (MID, UID))
        conn.execute("INSERT INTO candidate VALUES (?, 7, 'model', 'gain', 1, 0.5)", (MID,))
        conn.executemany("INSERT INTO grade VALUES (?, 7, ?, ?, 'no_card_hit', '', 2, "
                         "'prevalence;cpss', 'candidate')",
                         [(MID, "allele_aware", 1), (MID, "homolog_only", 0)])
        conn.execute("INSERT INTO release VALUES ('1.0.0-test', ?, '2026-11-03T09:00:00+00:00', "
                     "?, 'abc123', 0, '4.0.1', '{\"pyseer\": \"1.3.12\"}', ?, ?)",
                     (kb.SCHEMA_VERSION, "b" * 64, license, doi))
    conn.close()
    (path.parent / "kb_report.json").write_text(json.dumps({"sha256": sha256_file(path)}))
    return path


def _packaged(tmp_path: Path, name="kb", **release):
    d = tmp_path / name
    d.mkdir()
    db = _kb(d / "kanit.sqlite", **release)
    return db, m.package(db, d / "package")


def _zip_members(out: Path) -> dict[str, bytes]:
    with zipfile.ZipFile(out / m.TABLES_ZIP) as z:
        return {n: z.read(n) for n in z.namelist()}


def test_the_dictionary_describes_exactly_the_schema():
    d = _dictionary()
    assert d["schema_version"] == kb.SCHEMA_VERSION
    conn = sqlite3.connect(":memory:")
    kb.create(conn)
    tables, views = m.describe(conn, d)
    assert len(tables) == len(d["tables"]) and {v["name"] for v in views} == set(d["views"])
    grade = next(t for t in tables if t["name"] == "grade")
    cols = {c["name"]: c for c in grade["columns"]}
    assert cols["grade"]["enum"] == ["confirmed", "strong_novel", "candidate", "weak", "none"]
    assert cols["n_layers"]["range"] == (0, 4) and cols["primary_mode"]["enum"] == [0, 1]
    assert grade["foreign_keys"] == [{"fields": ["model_id", "pattern_id"],
                                      "resource": "candidate",
                                      "reference": ["model_id", "pattern_id"]}]
    aro = next(t for t in tables if t["name"] == "aro")    # a CHECK over two lines
    assert next(c for c in aro["columns"] if c["name"] == "model_type")["enum"] == [
        "homolog", "variant", "overexpression", "rrna_variant"]
    snap = next(t for t in tables if t["name"] == "data_snapshot")   # REFERENCES organism
    assert snap["foreign_keys"][0]["reference"] == ["organism_id"]


def test_a_dictionary_that_does_not_match_the_schema_stops_it():
    d = copy.deepcopy(_dictionary())
    del d["tables"]["grade"]["columns"]["layers_passed"]
    d["tables"]["grade"]["columns"]["stale"] = "a column no longer there"
    del d["views"]["v_novel"]
    conn = sqlite3.connect(":memory:")
    kb.create(conn)
    with pytest.raises(ValueError) as e:
        m.describe(conn, d)
    for part in ("grade.layers_passed is not described", "grade.stale is described",
                 "view v_novel is not described"):
        assert part in str(e.value)


def test_a_draft_package_holds_the_deposit_and_its_checksums(tmp_path):
    db, r = _packaged(tmp_path)
    out = Path(r["directory"])
    assert sorted(p.name for p in out.iterdir()) == sorted([*m.FILES, "SHA256SUMS"])
    for line in (out / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ")
        assert sha256_file(out / name) == digest
    assert sha256_file(out / "kanit.sqlite") == sha256_file(db)
    assert r["draft"] and r["missing"] == ["license", "doi"]
    assert r["tables"]["grade"] == 2 and r["tables"]["external_call"] == 0
    assert json.loads((db.parent / "package_report.json").read_text()) == r
    assert "Draft package" in (out / "README.md").read_text()
    assert "Draft package" in (out / "DATA_DICTIONARY.md").read_text()
    cff = yaml.safe_load((out / "CITATION.cff").read_text())
    assert cff["type"] == "dataset" and cff["version"] == "1.0.0-test"
    assert cff["date-released"] == "2026-11-03" and "doi" not in cff and "license" not in cff
    assert cff["authors"] and cff["title"].startswith("KANIT")


def test_the_tables_come_back_from_the_csv(tmp_path):
    db, r = _packaged(tmp_path)
    files = _zip_members(Path(r["directory"]))
    dp = json.loads(files["datapackage.json"])
    assert "DATA_DICTIONARY.md" in files and "licenses" not in dp and "id" not in dp
    conn = sqlite3.connect(db)
    assert {x["name"] for x in dp["resources"]} == {
        n for (n,) in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
    for res in dp["resources"]:
        raw = files[res["path"]]
        assert res["bytes"] == len(raw) and res["hash"] == "sha256:" + hashlib.sha256(raw).hexdigest()
        rows = list(csv.reader(io.StringIO(raw.decode("utf-8"))))
        assert rows[0] == [f["name"] for f in res["schema"]["fields"]]
        n = conn.execute(f"SELECT count(*) FROM {res['name']}").fetchone()[0]
        assert len(rows) - 1 == n == r["tables"][res["name"]]
    by = {x["name"]: x for x in dp["resources"]}

    def field(table, name):
        return next(f for f in by[table]["schema"]["fields"] if f["name"] == name)

    pattern = list(csv.DictReader(io.StringIO(files["tables/pattern.csv"].decode())))
    assert [p["pattern_id"] for p in pattern] == ["7", "8"]          # primary key order
    assert base64.b64decode(pattern[0]["carriers"]) == CARRIERS and pattern[1]["carriers"] == ""
    assert field("pattern", "carriers")["format"] == "binary"
    lp = next(csv.DictReader(io.StringIO(files["tables/label_permutation.csv"].decode())))
    assert lp["flag"] == "" and lp["null_sd"] == "" and lp["p"] == "0.0099"
    # NOT NULL text that holds the empty string is not required; its enum drops ''
    assert field("label_permutation", "flag")["constraints"] == {
        "enum": ["permutation_not_significant"]}
    assert "constraints" not in field("grade", "card_reasons")
    assert field("grade", "grade")["constraints"]["required"] is True
    assert field("grade", "n_layers")["constraints"] == {"required": True, "minimum": 0,
                                                         "maximum": 4}
    assert by["grade"]["schema"]["primaryKey"] == ["model_id", "pattern_id", "mode"]
    assert {"fields": ["model_id", "pattern_id"], "reference": {
        "resource": "candidate", "fields": ["model_id", "pattern_id"]}} in \
        by["grade"]["schema"]["foreignKeys"]


def test_a_release_package_carries_its_licence_and_doi(tmp_path):
    _, r = _packaged(tmp_path, license="CC-BY-4.0", doi="10.5281/zenodo.1")
    out = Path(r["directory"])
    assert not r["draft"] and r["missing"] == []
    cff = yaml.safe_load((out / "CITATION.cff").read_text())
    assert cff["doi"] == "10.5281/zenodo.1" and cff["license"] == "CC-BY-4.0"
    dp = json.loads(_zip_members(out)["datapackage.json"])
    assert dp["id"] == "https://doi.org/10.5281/zenodo.1"
    assert dp["licenses"][0]["path"] == "https://creativecommons.org/licenses/by/4.0/"
    text = (out / "README.md").read_text()
    assert "Draft" not in text and "CC-BY-4.0" in text and "https://doi.org/10.5281/zenodo.1" in text


def test_the_same_database_gives_the_same_package(tmp_path):
    _, a = _packaged(tmp_path, "a")
    d = tmp_path / "b"
    d.mkdir()
    (d / "kanit.sqlite").write_bytes((tmp_path / "a" / "kanit.sqlite").read_bytes())
    (d / "kb_report.json").write_text((tmp_path / "a" / "kb_report.json").read_text())
    b = m.package(d / "kanit.sqlite", d / "package")
    assert a["files"] == b["files"]


def test_a_database_the_build_did_not_check_is_refused(tmp_path):
    d = tmp_path / "kb"
    d.mkdir()
    db = _kb(d / "kanit.sqlite")
    (d / "kb_report.json").write_text(json.dumps({"sha256": "0" * 64}))
    with pytest.raises(ValueError, match="not the file the build checked"):
        m.package(db, d / "package")
    assert not (d / "package").exists() and not (d / "package.tmp").exists()


def test_the_package_validates_as_a_data_package(tmp_path):
    frictionless = pytest.importorskip("frictionless")
    _, r = _packaged(tmp_path)
    unpacked = tmp_path / "unpacked"
    with zipfile.ZipFile(Path(r["directory"]) / m.TABLES_ZIP) as z:
        z.extractall(unpacked)
    report = frictionless.validate(str(unpacked / "datapackage.json"))
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_the_readme_does_not_repeat_the_tool_name():
    assert m._bare_version("poppunk", "poppunk 2.7.8") == "2.7.8"
    assert m._bare_version("unitig_caller", "unitig-caller 1.3.2") == "1.3.2"
    assert m._bare_version("quast", "QUAST v5.3.0") == "v5.3.0"
    assert m._bare_version("checkm2", "1.1.0") == "1.1.0"
    assert m._bare_version("graph_tool", "3.0 (commit 1a3c1a29)") == "3.0 (commit 1a3c1a29)"

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Guard: every declared version must match lib.kb_schema.KB_SCHEMA_VERSION.

Five files carry the project version, and copies drift. .zenodo.json matters
most: it is the text Zenodo mints a permanent DOI over.

kb_api imports the constant so it cannot drift. Static files (.zenodo.json,
CITATION.cff, pyproject.toml, config.yaml) cannot import anything, so this test
is their equivalent: bump the schema and CI tells you what else to bump.
"""

import json
import re
import sys
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.kb_schema import KB_SCHEMA_VERSION  # noqa: E402


def _declared():
    root = PROJECT_ROOT
    return {
        ".zenodo.json": json.loads((root / ".zenodo.json").read_text())["version"],
        "CITATION.cff": yaml.safe_load((root / "CITATION.cff").read_text())["version"],
        "pyproject.toml": re.search(r'^version = "(.*)"',
                                    (root / "pyproject.toml").read_text(), re.M).group(1),
        "config.yaml": yaml.safe_load(
            (root / "config" / "config.yaml").read_text())["project"]["version"],
    }


def test_all_declared_versions_match_the_schema():
    drift = {f: v for f, v in _declared().items() if v != KB_SCHEMA_VERSION}
    assert not drift, (
        f"version drift vs KB_SCHEMA_VERSION={KB_SCHEMA_VERSION}: {drift}. "
        ".zenodo.json in particular gets a permanent DOI — fix before depositing.")


def test_no_stale_schema_version_in_public_prose():
    """Version fields are not the only place a version hides.

    .zenodo.json's `description` and `notes` can state the schema in running
    text, and the prose is what a human reads off the DOI landing page. Catch any
    schema-like number in that prose that is not the current one.
    """
    z = json.loads((PROJECT_ROOT / ".zenodo.json").read_text())
    prose = f"{z.get('description', '')} {z.get('notes', '')}"
    mentioned = set(re.findall(r"schema (\d+\.\d+\.\d+)", prose))
    stale = mentioned - {KB_SCHEMA_VERSION}
    assert not stale, (
        f".zenodo.json prose mentions schema {sorted(stale)} but the schema is "
        f"{KB_SCHEMA_VERSION}. This text gets a permanent DOI.")


def test_kb_schema_version_is_not_duplicated_in_config():
    """config.yaml must not carry its own copy of the schema version.

    A copy there drifts, and run_metadata.collect_versions would copy the stale
    value into the KB. The constant in lib/kb_schema is the one source.
    """
    cfg = yaml.safe_load((PROJECT_ROOT / "config" / "config.yaml").read_text())
    prov = cfg.get("provenance", {}) or {}
    assert "kb_schema_version" not in prov, (
        "config.yaml provenance.kb_schema_version is back — it will drift from "
        "lib.kb_schema.KB_SCHEMA_VERSION and poison run_metadata. Read the code "
        "constant instead.")


def test_collect_versions_reports_the_code_schema_version():
    from lib import run_metadata
    v = run_metadata.collect_versions(config={"provenance": {"kb_schema_version": "0.0.1"}})
    assert v["kb_schema_version"] == KB_SCHEMA_VERSION   # config cannot override it

#!/usr/bin/env python3
"""The data of the end-to-end test (scripts/e2e_subset.py) and the config overlay."""
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import config as cfgmod  # noqa: E402

pytestmark = pytest.mark.unit


def test_overlay_and_prefix(tmp_path, monkeypatch):
    over = tmp_path / "over.yaml"
    over.write_text("paths_prefix: e2e\npaths_prefix_exclude: [rgi_db_dir]\ncv: {n_folds: 3}\n")
    monkeypatch.setenv(cfgmod.OVERLAY_ENV, str(over))
    c = cfgmod.load_config()
    assert c["cv"]["n_folds"] == 3 and c["cv"]["n_repeats"] == 5      # merged, not replaced
    assert cfgmod.resolve_path("cv_dir", organism="o", antibiotic="a", config=c) == \
        PROJECT_ROOT / "e2e/results/o/a/cv"
    assert cfgmod.resolve_path("rgi_db_dir", config=c) == PROJECT_ROOT / "data/external/rgi_db"
    assert "paths_prefix" not in cfgmod.load_config(cfgmod.CONFIG_FILE)  # explicit path: no overlay


def test_subset(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("amrtest_e2e", PROJECT_ROOT / "scripts" /
                                                  "e2e_subset.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    monkeypatch.setattr(cfgmod, "PROJECT_ROOT", tmp_path)
    base = cfgmod.load_config(cfgmod.CONFIG_FILE)
    over = yaml.safe_load((PROJECT_ROOT / "tests/data/e2e_overlay.yaml").read_text())
    dst_cfg = cfgmod.deep_merge(base, over)
    meta = cfgmod.resolve_path("metadata_file", organism="ecoli", config=base).parent
    raw = cfgmod.resolve_path("raw_genomes_dir", organism="ecoli", config=base)
    meta.mkdir(parents=True)
    raw.mkdir(parents=True)
    ids = [f"562.{i}" for i in range(40)]
    long = pd.DataFrame({"genome_id": ids * 2, "antibiotic": ["ampicillin"] * 40 + ["colistin"] * 40,
                         "label": [i % 2 for i in range(40)] * 2})
    long.to_csv(meta / "amr_cleaned_long.csv", index=False)
    pd.DataFrame({"genome_id": ids, "status": ["passed"] * 39 + ["failed"]}).to_csv(
        meta / "download_report.csv", index=False)
    for name in ("genomes.csv", "amr_records.csv"):
        pd.DataFrame({"genome_id": ids, "x": 1}).to_csv(meta / name, index=False)
    (meta / "query.json").write_text("{}")
    (meta / "snapshot.json").write_text(json.dumps({"query": {"api_version": "1.9.3"}}))
    for g in ids:
        (raw / f"{g}.fna").write_text(">c\nACGT\n")
    picked = m.subset(base, dst_cfg, "ecoli", "ampicillin", 20, seed=1)
    assert len(picked) == 20 and "562.39" not in picked          # failed assembly left out
    dst = cfgmod.resolve_path("metadata_file", organism="ecoli", config=dst_cfg)
    assert str(dst).startswith(str(tmp_path / "e2e"))
    wide = pd.read_csv(dst, dtype={"Genome ID": str})
    assert list(wide.columns) == ["Genome ID", "ampicillin"] and wide["ampicillin"].sum() == 10
    link = cfgmod.resolve_path("raw_genomes_dir", organism="ecoli", config=dst_cfg) / f"{picked[0]}.fna"
    assert link.is_symlink() and link.read_text() == ">c\nACGT\n"
    snap = json.loads((dst.parent / "snapshot.json").read_text())
    assert snap["e2e_subset"]["antibiotic"] == "ampicillin" and snap["n_genomes"] == 20

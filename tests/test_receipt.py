#!/usr/bin/env python3
"""The receipt of a workflow task, which is also its step manifest (scripts/receipt.py)."""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
pytestmark = pytest.mark.unit


def test_receipt_records_the_task(tmp_path):
    (tmp_path / ".command.begin").touch()
    os.utime(tmp_path / ".command.begin", (1_790_000_000, 1_790_000_000))
    dep = tmp_path / "dep1.json"
    dep.write_text(json.dumps({"step": "cv_folds", "key": "ecoli__ampicillin"}))
    env = {**os.environ, "SLURM_JOB_ID": "6443221", "SLURMD_NODENAME": "hamsi20",
           "APPTAINER_CONTAINER": "/arf/scratch/u/amr/containers/amr.sif"}
    out = subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "scripts" / "receipt.py"), "--step", "cv_unit",
         "--key", "ecoli__ampicillin lineage_aware 1 0", "--process", "MODELS:CV_UNIT",
         "--attempt", "2", "--cpus", "20", "--memory", "128 GB",
         "--tool", "echo=echo tool 1.2.3", "--tool", "missing=no_such_tool_kanit --version"],
        cwd=tmp_path, env=env, capture_output=True, text=True, check=True).stdout
    r = json.loads(out)
    assert (r["step"], r["process"], r["attempt"], r["cpus"], r["memory"]) == (
        "cv_unit", "MODELS:CV_UNIT", 2, 20, "128 GB")
    assert r["started_at"] == "2026-09-21T14:13:20+00:00" and r["seconds"] > 0
    assert r["slurm"] == {"job_id": "6443221", "node": "hamsi20"}
    assert r["environment"]["container"].endswith("amr.sif") and r["environment"]["python"]
    assert r["inputs"] == [{"file": "dep1.json", "step": "cv_folds", "key": "ecoli__ampicillin",
                            "sha256": hashlib.sha256(dep.read_bytes()).hexdigest()}]
    assert r["tools"] == {"echo": "tool 1.2.3", "missing": None}
    assert r["work_dir"] == str(tmp_path)


def test_receipt_without_nextflow(tmp_path):
    out = subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "receipt.py"),
                          "--step", "panel", "--stub"], cwd=tmp_path, capture_output=True,
                         text=True, check=True).stdout
    r = json.loads(out)
    assert r["stub"] and r["started_at"] is None and r["inputs"] == [] and r["tools"] == {}

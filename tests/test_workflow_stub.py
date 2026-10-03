#!/usr/bin/env python3
"""The Nextflow skeleton connects every step: -stub-run of each entry (needs
Nextflow; runs when KANIT_RUN_NEXTFLOW=1, since it takes about a minute)."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
NEXTFLOW = os.environ.get("NEXTFLOW", shutil.which("nextflow") or "")

pytestmark = pytest.mark.skipif(
    os.environ.get("KANIT_RUN_NEXTFLOW") != "1" or not NEXTFLOW,
    reason="set KANIT_RUN_NEXTFLOW=1 (and NEXTFLOW if nextflow is not on PATH)")


def _run(tmp_path, *args):
    cmd = [NEXTFLOW, "-q", "run", str(PROJECT_ROOT / "main.nf"), "-stub-run", "-profile", "test",
           "--python", sys.executable, "--trace_dir", str(tmp_path / "trace"), *args]
    r = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
    return r, tmp_path / "trace" / "trace.tsv"


def test_default_entry_reaches_every_step(tmp_path):
    r, trace = _run(tmp_path, "--organisms", "ecoli,kpneumoniae", "--stub_panel",
                    str(PROJECT_ROOT / "tests" / "data" / "stub_panel_decisions.csv"))
    assert r.returncode == 0, r.stdout + r.stderr
    t = pd.read_csv(trace, sep="\t")
    n = t["process"].str.split(":").str[-1].value_counts().to_dict()
    assert (t["status"] == "COMPLETED").all()
    assert n["CV_UNIT"] == 3 * 50 and n["LP_CHUNK"] == 3 * 50 and n["CPSS_CHUNK"] == 3 * 10
    for p in ("PANEL", "GRADING", "LP_ACROSS", "EXTERNAL_COMPARE", "CARD_LAYER", "PYSEER_POST"):
        assert n.get(p), p


def test_parameters_are_checked(tmp_path):
    r, _ = _run(tmp_path)
    assert r.returncode != 0 and "organisms" in (r.stdout + r.stderr)
    r, _ = _run(tmp_path, "--organisms", "salmonella")
    assert r.returncode != 0 and "Unknown organism" in (r.stdout + r.stderr)

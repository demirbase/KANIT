#!/usr/bin/env python3
"""The Nextflow skeleton connects every step: -stub-run of each entry (needs
Nextflow; runs when KANIT_RUN_NEXTFLOW=1, since it takes about a minute)."""
import json
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
    # 50 outer folds and 50 permutation chunks per model, 4 per job in the stub (lib/packing.py)
    assert n["CV_UNITS"] == 3 * 13 and n["LP_CHUNKS"] == 3 * 13 and n["CPSS_CHUNK"] == 3 * 1
    for p in ("PANEL", "GRADING", "LP_ACROSS", "EXTERNAL_COMPARE", "CARD_LAYER", "PYSEER_POST"):
        assert n.get(p), p
    # the patterns nominated by association (protocol §14 item 6): every layer and the
    # grades of every model again
    for p in ("ASSOCIATION_SET", "PREVALENCE_ASSOCIATION", "MDA_ASSOCIATION",
              "CARD_LAYER_ASSOCIATION", "GRADING_ASSOCIATION"):
        assert n.get(p) == n["GRADING"], p
    assert n["RUN_OUTPUTS"] == n["SACCT_DUMP"] == n["RUN_RESOURCES"] == 1
    assert n["REPORTS"] == n["COMPLETENESS"] == n["COMPLETENESS_GATE"] == 1
    assert "BACKUP_PACK" not in n                     # no --backup_remote: no backup
    # the unitig store reads the panel: every store starts after the panel is written
    proc = t["process"].str.split(":").str[-1]
    assert t.loc[proc == "STORE", "start"].min() >= t.loc[proc == "PANEL", "complete"].max()
    # the run manifest and every task's receipt (its step manifest) are kept with the run
    run = json.loads((tmp_path / "trace" / "run_manifest.json").read_text())
    assert run["run"]["entry"] == "main" and run["run"]["stub"] is True
    assert len(run["code"]["commit"]) == 40 and run["completed"]["success"] is True
    assert run["config"]["config_yaml"]["sha256"] and run["config"]["protocol"]["version"] == "1.0"
    assert set(run["snapshots"]) == {"ecoli", "kpneumoniae"}
    receipts = list((tmp_path / "trace" / "tasks").glob("*/*/receipt.json"))
    assert len(receipts) == len(t)
    r = json.loads((tmp_path / "trace" / "tasks" / "PANEL" / "run" / "receipt.json").read_text())
    assert r["step"] == "panel" and r["process"].endswith("PANEL") and r["stub"] is True
    keys = {json.loads(f.read_text())["key"]
            for f in (tmp_path / "trace" / "tasks").glob("GRADING*/*/receipt.json")}
    assert len(keys) == 2 * n["GRADING"]
    assert sum(k.endswith(" association") for k in keys) == n["GRADING"]


def test_parameters_are_checked(tmp_path):
    r, _ = _run(tmp_path)
    assert r.returncode != 0 and "organisms" in (r.stdout + r.stderr)
    r, _ = _run(tmp_path, "--organisms", "salmonella")
    assert r.returncode != 0 and "Unknown organism" in (r.stdout + r.stderr)


@pytest.mark.parametrize("entry, expected", [
    ("DOWNLOAD", {"DOWNLOAD_BVBRC": 2}),
    ("CONTEXT", {"CONTEXT_QUERY": 2, "CONTEXT_BUILD": 2})])
def test_internet_entries(tmp_path, entry, expected):
    expected = {**expected, "REPORTS": 1, "COMPLETENESS": 1, "RUN_OUTPUTS": 1, "BACKUP_PACK": 1,
                "BACKUP_UPLOAD": 1, "SACCT_DUMP": 1, "RUN_RESOURCES": 1, "COMPLETENESS_GATE": 1}
    r, trace = _run(tmp_path, "-entry", entry, "--organisms", "ecoli,kpneumoniae",
                    "--backup_remote", str(tmp_path / "remote"))
    assert r.returncode == 0, r.stdout + r.stderr
    t = pd.read_csv(trace, sep="\t")
    assert (t["status"] == "COMPLETED").all()
    assert t["process"].str.split(":").str[-1].value_counts().to_dict() == expected


def test_cabbage_entry(tmp_path):
    """External validation: the global steps once, every per-organism step once per
    organism, the tools once per shard."""
    r, trace = _run(tmp_path, "-entry", "CABBAGE", "--organisms", "ecoli,kpneumoniae",
                    "--rgi_shards", "2", "--external_shards", "3")
    assert r.returncode == 0, r.stdout + r.stderr
    t = pd.read_csv(trace, sep="\t")
    assert (t["status"] == "COMPLETED").all()
    n = t["process"].str.split(":").str[-1].value_counts().to_dict()
    assert {k: n.get(k) for k in ("CABBAGE_DOWNLOAD", "CABBAGE_SELECT", "CABBAGE_FETCH")} == \
        {"CABBAGE_DOWNLOAD": 1, "CABBAGE_SELECT": 1, "CABBAGE_FETCH": 1}
    for p in ("QC_PATHS", "CHECKM2", "PREPARE", "ASSIGN", "CALL", "PREDICT", "EXTERNAL_PREP",
              "EXTERNAL_COLLECT", "RGI_COLLECT", "COMPARE"):
        assert n.get(f"CABBAGE_{p}") == 2, p
    assert n["CABBAGE_EXTERNAL_RUN"] == 6 and n["CABBAGE_RGI"] == 4
    assert n["COMPLETENESS_GATE"] == 1


def test_internet_tasks_run_on_the_login_node(tmp_path):
    """TRUBA's compute nodes have no internet: the tasks that need it run where the
    head job runs."""
    r = subprocess.run([NEXTFLOW, "config", "-flat", "-profile", "truba", str(PROJECT_ROOT)],
                       cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "process.'withLabel:internet'.executor = 'local'" in r.stdout
    for f in ("download.nf", "context.nf"):
        text = (PROJECT_ROOT / "modules" / "local" / f).read_text()
        assert text.count("label 'internet'") == {"download.nf": 1, "context.nf": 1}[f]


def test_python_output_is_unbuffered_on_truba(tmp_path):
    """The logs of a running job show its progress (the units of a packed job wrote nothing
    until they ended). Set in the env scope, which leaves the tasks' cache keys unchanged."""
    r = subprocess.run([NEXTFLOW, "config", "-flat", "-profile", "truba", str(PROJECT_ROOT)],
                       cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "env.PYTHONUNBUFFERED = '1'" in r.stdout

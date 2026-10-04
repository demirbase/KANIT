#!/usr/bin/env python3
"""Resources of every task of a run (scripts/run_resources.py)."""
import importlib.util
import json
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
pytestmark = pytest.mark.unit


def _script():
    spec = importlib.util.spec_from_file_location("amrtest_run_resources",
                                                  PROJECT_ROOT / "scripts" / "run_resources.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_parsers():
    m = _script()
    assert m.seconds("1-02:03:04") == 93784 and m.seconds("05:30.5") == 330.5
    assert m.seconds("") is None and m.size_bytes("1.5G") == 1.5 * 1024 ** 3
    assert m.size_bytes("4000Mc") == 4000 * 1024 ** 2 and m.size_bytes("0") == 0


def test_receipts_trace_and_sacct_are_joined(tmp_path):
    m = _script()
    for proc, tag, job in (("CV_UNIT", "a", "101"), ("CV_UNIT", "b", "102"), ("PANEL", "run", None)):
        d = tmp_path / "tasks" / proc / tag
        d.mkdir(parents=True)
        (d / "receipt.json").write_text(json.dumps({
            "process": f"MODELS:{proc}", "step": proc.lower(), "key": tag, "attempt": 1,
            "cpus": 20, "slurm": {"job_id": job, "node": "barbun1"}}))
    (tmp_path / "sacct.tsv").write_text(
        "JobID|JobName|Partition|State|ExitCode|Elapsed|Start|End|AllocCPUS|ReqMem|MaxRSS|TotalCPU|NodeList\n"
        "101|kanit_CV_UNIT_a|barbun|COMPLETED|0:0|02:00:00|s|e|20|190G||30:00:00|barbun1\n"
        "101.batch|batch||COMPLETED|0:0|02:00:00|s|e|20||80G|30:00:00|barbun1\n"
        "101.extern|extern||COMPLETED|0:0|02:00:00|s|e|20||1M|00:00:01|barbun1\n"
        "102|kanit_CV_UNIT_b|barbun|COMPLETED|0:0|01:00:00|s|e|20|190G||10:00:00|barbun2\n"
        "102.batch|batch||COMPLETED|0:0|01:00:00|s|e|20||100G|10:00:00|barbun2\n")
    (tmp_path / "trace.tsv").write_text("task_id\tnative_id\tprocess\trealtime\t%cpu\tpeak_rss\n"
                                        "1\t101\tMODELS:CV_UNIT\t7200000\t1490.5\t70000000000\n")
    table, summary = m.resources(tmp_path)
    cv = summary["by_process"]["MODELS:CV_UNIT"]
    assert (cv["tasks"], cv["core_hours"], cv["cpu_hours"], cv["max_rss_gib"]) == (2, 60.0, 40.0,
                                                                                     100.0)
    assert cv["max_trace_peak_rss_gib"] == round(70e9 / 1024 ** 3, 2)
    assert summary["by_process"]["MODELS:PANEL"]["core_hours"] == 0
    assert summary["jobs_without_accounting"] == 0 and len(table) == 3

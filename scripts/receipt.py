#!/usr/bin/env python3
"""Receipt of one workflow task (stdout, JSON), which is also the task's step manifest.

The steps write their outputs into the shared results tree; a receipt is what a
Nextflow task hands to the tasks that depend on it, so that a rerun upstream reruns
them too. Every receipt is also published to the run directory (nextflow.config).

It records the step and its key; the process, attempt and resources Nextflow gave
the task; the start (the task's .command.begin) and the end; the host, SLURM job and
node; the code commit; the checksum, step and key of every receipt the task took as
input (dep*.json); the container, Python and the versions of the analysis packages;
and the version of every tool named with --tool. CPU time and peak memory are in
Nextflow's trace. Standard library only: it also runs in the tool containers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

try:
    from importlib import metadata
except ImportError:                       # Python 3.7 in a tool container
    metadata = None  # type: ignore[assignment]

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PACKAGES = ["numpy", "pandas", "scipy", "scikit-learn", "xgboost", "statsmodels", "PyYAML",
            "poppunk", "unitig-caller", "pyseer", "biopython", "rgi", "CheckM2", "quast",
            "resfinder"]


def _commit():
    try:
        return subprocess.run(["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _iso(ts: float) -> str:
    return datetime.fromtimestamp(ts, timezone.utc).isoformat(timespec="seconds")


def _inputs(cwd: Path) -> list:
    out = []
    for f in sorted(cwd.glob("dep*.json")):
        data = f.read_bytes()
        try:
            r = json.loads(data)
        except ValueError:
            r = {}
        out.append({"file": f.name, "sha256": hashlib.sha256(data).hexdigest(),
                    "step": r.get("step"), "key": r.get("key")})
    return out


def _packages() -> dict:
    out: dict = {}
    if metadata is None:
        return out
    for p in PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            pass
    return out


def tool_version(command: str) -> str | None:
    """First non-empty line a version command prints (stdout, else stderr); None when
    the command fails."""
    try:
        r = subprocess.run(command, shell=True, capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if r.returncode != 0:
        return None
    for text in (r.stdout, r.stderr):
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        if lines:
            return lines[0]
    return None


def receipt(args, cwd: Path, env=os.environ) -> dict:
    now = datetime.now(timezone.utc)
    begin = cwd / ".command.begin"
    started = begin.stat().st_mtime if begin.exists() else None
    tools = {}
    for spec in args.tool:
        name, _, command = spec.partition("=")
        tools[name] = tool_version(command)
    return {
        "step": args.step, "key": args.key, "stub": args.stub,
        "process": args.process, "attempt": args.attempt, "cpus": args.cpus,
        "memory": args.memory or None,
        "started_at": _iso(started) if started else None,
        "finished_at": now.isoformat(timespec="seconds"),
        "seconds": round(now.timestamp() - started, 1) if started else None,
        "host": platform.node(),
        "slurm": {"job_id": env.get("SLURM_JOB_ID"),
                  "node": env.get("SLURMD_NODENAME") or env.get("SLURM_JOB_NODELIST")},
        "work_dir": str(cwd),
        "code_commit": _commit(),
        "inputs": _inputs(cwd),
        "environment": {"container": env.get("APPTAINER_CONTAINER")
                        or env.get("SINGULARITY_CONTAINER"),
                        "python": platform.python_version(), "packages": _packages()},
        "tools": tools,
    }


def main():
    ap = argparse.ArgumentParser(description="Receipt (step manifest) of one workflow task.")
    ap.add_argument("--step", required=True)
    ap.add_argument("--key", default="")
    ap.add_argument("--stub", action="store_true", help="written by a -stub-run task")
    ap.add_argument("--process", default=None, help="Nextflow process name")
    ap.add_argument("--attempt", type=int, default=None)
    ap.add_argument("--cpus", type=int, default=None)
    ap.add_argument("--memory", default=None)
    ap.add_argument("--tool", action="append", default=[], metavar="NAME=COMMAND",
                    help="record the version a command prints (repeatable)")
    args = ap.parse_args()
    json.dump(receipt(args, Path.cwd()), sys.stdout, sort_keys=True)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()

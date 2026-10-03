#!/usr/bin/env python3
"""Receipt of one workflow task (stdout, JSON): the step, its key, the time, the host
and the code commit. The steps write their outputs into the shared results tree;
a receipt is what a Nextflow task hands to the tasks that depend on it, so that a
rerun upstream reruns them too. Standard library only: it also runs in the tool
containers."""
import argparse
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _commit():
    try:
        return subprocess.run(["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main():
    ap = argparse.ArgumentParser(description="Receipt of one workflow task.")
    ap.add_argument("--step", required=True)
    ap.add_argument("--key", default="")
    ap.add_argument("--stub", action="store_true", help="written by a -stub-run task")
    args = ap.parse_args()
    print(json.dumps({"step": args.step, "key": args.key, "stub": args.stub,
                      "finished_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                      "host": platform.node(), "code_commit": _commit()}, sort_keys=True))


if __name__ == "__main__":
    main()

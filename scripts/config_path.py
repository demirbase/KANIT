#!/usr/bin/env python3
"""Print one path of config.yaml's paths_organism block, resolved (for workflow tasks)."""
import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config, resolve_path  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description="Resolve a config path.")
    ap.add_argument("key")
    ap.add_argument("--organism")
    ap.add_argument("--antibiotic")
    args = ap.parse_args()
    print(resolve_path(args.key, organism=args.organism, antibiotic=args.antibiotic,
                       config=load_config()))


if __name__ == "__main__":
    main()

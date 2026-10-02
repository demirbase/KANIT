#!/usr/bin/env python3
"""Step 10 — prevalence layer of one panel pair's candidate patterns (protocol §8.2).

Prevalence among the resistant and the susceptible genomes of the model, Δ and
its direction, and a Fisher exact test with Benjamini–Hochberg over the
candidates. The rule is in lib/prevalence.py.

Inputs: candidates_file, the model matrix. Output: layers_dir/prevalence.csv.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import get_target, load_config, resolve_path  # noqa: E402
from lib.matrix_store import ModelMatrix  # noqa: E402
from lib.prevalence import prevalence_layer  # noqa: E402


def main():
    config = load_config()
    default_org, default_ab = get_target(config=config)
    ap = argparse.ArgumentParser(description="Prevalence layer of one panel pair's candidates.")
    ap.add_argument("--organism", default=default_org)
    ap.add_argument("--antibiotic", default=default_ab)
    args = ap.parse_args()
    org, ab, cfg = args.organism, args.antibiotic, config["prevalence"]

    def path(key):
        return resolve_path(key, organism=org, antibiotic=ab, config=config)

    patterns = sorted(set(pd.read_csv(path("candidates_file"))["pattern_id"].astype(int)))
    layer = prevalence_layer(ModelMatrix(path("matrix_dir")), patterns,
                             min_delta=cfg["min_delta"], alpha=cfg["alpha"])
    out = path("layers_dir")
    out.mkdir(parents=True, exist_ok=True)
    layer.to_csv(out / "prevalence.csv", index=False)
    print(f"Prevalence — {org} / {ab}: {int(layer['passes'].sum())} of {len(layer)} "
          f"candidate patterns pass")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Peak memory of the in-memory XGBoost matrix (the storage rule of the protocol).

Builds the pattern matrix of every genome in an organism store, then a
QuantileDMatrix (max_bin = 2) from it in row batches, and reports the peak
resident memory. All genomes of the organism bound every model of that organism
from above; the labels are placeholders and do not change the matrix's size.

The rule: a peak of at most 300 GB means one file per model, otherwise row shards.

    python scripts/xgb_memory_probe.py --store-dir data/processed/ecoli/unitig_store \\
        --work-dir /path/to/scratch/probe_ecoli --out probe_ecoli.json
"""
import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import matrix_store  # noqa: E402
from lib.config import load_config  # noqa: E402
from lib.run_metadata import peak_rss_gb  # noqa: E402


class RowBatches(xgb.DataIter):
    """All rows of a model matrix, ``batch_rows`` at a time."""

    def __init__(self, mm, batch_rows):
        self.mm, self.batch_rows, self.i = mm, batch_rows, 0
        self.starts = list(range(0, mm.n_genomes, batch_rows))
        self.labels = mm.labels
        super().__init__()

    def next(self, input_data):
        if self.i == len(self.starts):
            return False
        s = self.starts[self.i]
        e = min(s + self.batch_rows, self.mm.n_genomes)
        input_data(data=self.mm.rows(np.arange(s, e)), label=self.labels[s:e])
        self.i += 1
        return True

    def reset(self):
        self.i = 0


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="Measure the XGBoost matrix's peak memory.")
    ap.add_argument("--store-dir", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True,
                    help="where the probe's pattern matrix is written")
    ap.add_argument("--out", type=Path, required=True, help="JSON report")
    ap.add_argument("--min-support", type=int, default=int(config["unitig"]["min_support"]))
    ap.add_argument("--batch-rows", type=int, default=128)
    ap.add_argument("--limit-gb", type=float, default=300.0)
    args = ap.parse_args()

    store = matrix_store.Store(args.store_dir)
    n = len(store.genomes)
    genomes = pd.DataFrame({"Genome ID": store.genomes, "label": np.arange(n) % 2,
                            "lineage": 0})
    t0 = time.time()
    s = matrix_store.build_model_matrix(args.store_dir, genomes, args.work_dir,
                                        min_support=args.min_support)
    t_matrix, peak_matrix = time.time() - t0, peak_rss_gb()
    print(f"pattern matrix: {s['n_genomes']:,} genomes × {s['n_patterns']:,} patterns "
          f"in {t_matrix:.0f} s, peak {peak_matrix:.1f} GB")

    mm = matrix_store.ModelMatrix(args.work_dir)
    t0 = time.time()
    dm = xgb.QuantileDMatrix(RowBatches(mm, args.batch_rows), max_bin=2)
    t_dm, peak = time.time() - t0, peak_rss_gb()
    print(f"QuantileDMatrix: {dm.num_row():,} × {dm.num_col():,} in {t_dm:.0f} s, "
          f"peak {peak:.1f} GB")

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "host": platform.node(), "xgboost_version": xgb.__version__,
        "store_dir": str(args.store_dir), "min_support": args.min_support,
        "n_genomes": s["n_genomes"], "n_unitigs_kept": s["n_unitigs_kept"],
        "n_patterns": s["n_patterns"], "batch_rows": args.batch_rows,
        "seconds_pattern_matrix": round(t_matrix, 1), "peak_gb_pattern_matrix": round(peak_matrix, 2),
        "seconds_quantile_dmatrix": round(t_dm, 1), "peak_gb": round(peak, 2),
        "limit_gb": args.limit_gb,
        "storage": "single file" if peak <= args.limit_gb else "row shards",
    }
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(f"storage rule: {report['storage']} (peak {peak:.1f} GB, limit {args.limit_gb:.0f} GB)")


if __name__ == "__main__":
    main()

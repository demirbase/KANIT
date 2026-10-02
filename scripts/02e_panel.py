#!/usr/bin/env python3
"""Step 02e — the panel.

For every organism in the registry, counts resistant and susceptible genomes per
antibiotic among the genomes that pass quality control (02d) and have a lineage
(02c), applies the panel rule (lib.panel) and writes one row per pair with its
decision and reason:

    results/panel/panel_decisions.csv
    results/panel/panel_summary.json

Every input must exist, and every antibiotic column must carry its canonical
registry name; otherwise the step stops instead of shrinking the panel. The step
exits with status 2 when a pair is blocked: a drug that passes
the rule but has no registry class has to be added to the antibiotic registry
before any model is trained.

Usage:
    python scripts/02e_panel.py                    # all registry organisms
    python scripts/02e_panel.py --organisms ecoli kpneumoniae
"""
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import panel, registry  # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402


def _inputs(organism, config):
    qc_dir = resolve_path("genome_qc_dir", organism=organism, config=config)
    lineage_dir = resolve_path("lineage_dir", organism=organism, config=config)
    return {
        "phenotypes": resolve_path("metadata_file", organism=organism, config=config),
        "qc_table": qc_dir / f"02d_genome_qc_{organism}.csv",
        "clusters": lineage_dir / "poppunk_clusters.csv",
    }


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="Apply the panel rule to every organism.")
    ap.add_argument("--organisms", nargs="+", default=list(registry.load_organisms()),
                    help="registry slugs (default: every organism in the registry)")
    ap.add_argument("--min-minority", type=int,
                    default=int(config["panel"]["min_minority"]),
                    help="smallest minority class that enters the panel (config panel.min_minority)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="output directory (default: config paths_organism.panel_dir)")
    args = ap.parse_args()

    unknown = [o for o in args.organisms if o not in registry.load_organisms()]
    if unknown:
        sys.exit(f"ERROR: not in the organism registry: {', '.join(unknown)}")

    inputs = {org: _inputs(org, config) for org in args.organisms}
    missing = [f"{org}: {name} {path}" for org, files in inputs.items()
               for name, path in files.items() if not Path(path).exists()]
    if missing:
        sys.exit("ERROR: missing inputs (run 00, 02c and 02d first):\n  " + "\n  ".join(missing))

    stale = []
    for org, files in inputs.items():
        header = pd.read_csv(files["phenotypes"], nrows=0, encoding="utf-8")
        stale += [f"{org}: '{name}' -> '{canonical}'"
                  for name, canonical in panel.non_canonical_columns(header)]
    if stale:
        sys.exit("ERROR: phenotype columns without their canonical registry name; "
                 "run 00 again with the current registry:\n  " + "\n  ".join(stale))

    print("=" * 78)
    print(f"PANEL — minority class >= {args.min_minority} after QC and lineage assignment")
    print("=" * 78)
    rows, organisms = [], {}
    for org, files in inputs.items():
        phenotypes = pd.read_csv(files["phenotypes"], dtype={"Genome ID": str}, encoding="utf-8")
        qc = pd.read_csv(files["qc_table"], dtype={"genome_id": str}, encoding="utf-8")
        clusters = pd.read_csv(files["clusters"], dtype={"Genome ID": str}, encoding="utf-8")
        eligible = panel.eligible_genomes(qc, clusters)
        org_rows = panel.pair_rows(org, phenotypes, eligible, min_minority=args.min_minority)
        rows.extend(org_rows)
        n_inc = sum(r["decision"] == panel.INCLUDED for r in org_rows)
        organisms[org] = {
            "n_phenotyped": int(len(phenotypes)),
            "n_qc_assessed": int(len(qc)),
            "n_qc_pass": int(qc["pass_overall"].astype(str).str.lower().eq("true").sum()),
            "n_with_lineage": int(clusters["Genome ID"].nunique()),
            "n_eligible": len(eligible),
            "n_eligible_phenotyped": int(phenotypes["Genome ID"].isin(eligible).sum()),
            "n_pairs": len(org_rows),
            "n_included": n_inc,
            "inputs": {name: {"path": str(Path(p).relative_to(PROJECT_ROOT))
                              if Path(p).is_relative_to(PROJECT_ROOT) else str(p),
                              "sha256": _sha256(p)} for name, p in files.items()},
        }
        print(f"  {org:<26} eligible {len(eligible):>6,} genomes | "
              f"{len(org_rows):>3} antibiotics | {n_inc:>3} included")

    df = pd.DataFrame(rows, columns=panel.COLUMNS)
    out_dir = args.out_dir or resolve_path("panel_dir", config=config)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "panel_decisions.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8")

    counts = df["decision"].value_counts().to_dict()
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "min_minority": args.min_minority,
        "n_pairs": int(len(df)),
        "n_included": int(counts.get(panel.INCLUDED, 0)),
        "n_excluded": int(counts.get(panel.EXCLUDED, 0)),
        "n_blocked": int(counts.get(panel.BLOCKED, 0)),
        "organisms": organisms,
    }
    (out_dir / "panel_summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                                encoding="utf-8")
    print("-" * 78)
    print(f"  {summary['n_included']} of {summary['n_pairs']} pairs included -> {csv_path}")

    blocked = df[df["decision"] == panel.BLOCKED]
    if not blocked.empty:
        print("\nBLOCKED — add these drugs to config/registry/antibiotics.yaml first:")
        for r in blocked.itertuples():
            print(f"  {r.organism} / {r.antibiotic} (minority {r.minority})")
        sys.exit(2)


if __name__ == "__main__":
    main()

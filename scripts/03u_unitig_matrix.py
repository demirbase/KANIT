#!/usr/bin/env python3
"""Step 03u — unitig store and model matrices.

Two modes:

  --build-store  ORGANISM level. Runs unitig-caller once over every genome of the
                 organism that enters at least one panel pair (k from config
                 unitig.k), and parses its presence table into the organism
                 store (lib.matrix_store.build_store).
  default        MODEL level, for one panel pair: the pair's genomes (pass QC,
                 have a lineage, tested for the antibiotic), the frequency filter
                 and the pattern collapse (lib.matrix_store.build_model_matrix).

Outputs go to config paths_organism.unitig_store_dir and matrix_dir; the formats
are described in lib/matrix_store.py. unitig-caller's output is logged to
<store>/call/unitig_caller.log.

Usage:
    python scripts/03u_unitig_matrix.py --organism ecoli --build-store
    python scripts/03u_unitig_matrix.py --organism ecoli --antibiotic ampicillin
"""
import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import matrix_store, panel  # noqa: E402
from lib.config import get_target, load_config, resolve_path, resolve_tool  # noqa: E402
from lib.io_utils import run_logged  # noqa: E402
from lib.run_metadata import _tool_version, write_versions  # noqa: E402


def _included(organism, config):
    panel_csv = resolve_path("panel_dir", config=config) / "panel_decisions.csv"
    if not panel_csv.exists():
        sys.exit(f"ERROR: no panel at {panel_csv}; run 02e first.")
    return [ab for org, ab in panel.included_pairs(panel_csv) if org == organism]


def store_genome_ids(phenotypes, qc_table, clusters, antibiotics):
    """Every genome of at least one of the pairs, sorted."""
    ids = set()
    for ab in antibiotics:
        ids |= set(panel.pair_genomes(phenotypes, qc_table, clusters, ab)["Genome ID"])
    return sorted(ids)


def call_unitigs(genome_ids, genomes_dir, work_dir, *, threads, k, tool):
    """unitig-caller over the genomes -> <work_dir>/unitigs.rtab.

    A finished call leaves unitigs.rtab.done; an Rtab without it is the remnant
    of an interrupted run and is never reused.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    rtab, done = work_dir / "unitigs.rtab", work_dir / "unitigs.rtab.done"
    if rtab.exists() and done.exists():
        print(f"  ✓ reusing the finished call: {rtab}")
        return rtab
    fna = [Path(genomes_dir) / f"{g}.fna" for g in genome_ids]
    missing = [str(p) for p in fna if not p.exists()]
    if missing:
        sys.exit(f"ERROR: {len(missing)} assembly file(s) missing, e.g. {missing[:3]}")
    refs = work_dir / "unitig_refs.txt"
    refs.write_text("".join(f"{p.resolve()}\n" for p in fna), encoding="utf-8")
    done.unlink(missing_ok=True)
    print(f"  Running unitig-caller on {len(fna):,} genomes (k={k}, threads={threads})...")
    run_logged([tool, "--call", "--refs", refs, "--rtab", "--kmer", k,
                "--out", work_dir / "unitigs", "--threads", threads],
               work_dir / "unitig_caller.log")
    if not rtab.exists():
        sys.exit(f"ERROR: unitig-caller did not write {rtab}")
    write_versions(work_dir / "versions.json", {"unitig_caller": _tool_version(tool)})
    done.write_text(datetime.now(timezone.utc).isoformat(timespec="seconds") + "\n")
    return rtab


def build_store(organism, config, *, threads, min_support, rtab=None):
    store_dir = resolve_path("unitig_store_dir", organism=organism, config=config)
    if rtab is None:
        antibiotics = _included(organism, config)
        if not antibiotics:
            sys.exit(f"ERROR: no {organism} pair is in the panel.")
        genome_ids = store_genome_ids(*panel.read_inputs(organism, config), antibiotics)
        tool = resolve_tool("unitig-caller")
        if not tool:
            sys.exit("ERROR: unitig-caller not found on PATH (or set AMR_UNITIG_CALLER_BIN).")
        rtab = call_unitigs(genome_ids, resolve_path("raw_genomes_dir", organism=organism,
                                                     config=config),
                            store_dir / "call", threads=threads,
                            k=int(config["unitig"]["k"]), tool=tool)
    print(f"  Parsing {rtab} ...")
    s = matrix_store.build_store(rtab, store_dir, min_support=min_support)
    print(f"  ✓ store {store_dir}: {s['n_unitigs']:,} of {s['n_unitigs_seen']:,} unitigs "
          f"× {s['n_genomes']:,} genomes (min_support={min_support})")


def build_model(organism, antibiotic, config, *, min_support):
    if antibiotic not in _included(organism, config):
        sys.exit(f"ERROR: {organism}/{antibiotic} is not in the panel (panel_decisions.csv).")
    genomes = panel.pair_genomes(*panel.read_inputs(organism, config), antibiotic)
    store_dir = resolve_path("unitig_store_dir", organism=organism, config=config)
    out_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic, config=config)
    s = matrix_store.build_model_matrix(store_dir, genomes, out_dir, min_support=min_support)
    print(f"  ✓ {out_dir}: {s['n_genomes']:,} genomes (R {s['n_resistant']:,} / "
          f"S {s['n_susceptible']:,}) × {s['n_patterns']:,} patterns from "
          f"{s['n_unitigs_kept']:,} unitigs")


def main():
    config = load_config()
    unitig_cfg = config.get("unitig", {}) or {}
    default_org, default_ab = get_target(config=config)
    ap = argparse.ArgumentParser(description="Build the unitig store or a model matrix.")
    ap.add_argument("--organism", default=default_org)
    ap.add_argument("--antibiotic", default=default_ab)
    ap.add_argument("--build-store", action="store_true",
                    help="organism level: run unitig-caller and build the store")
    ap.add_argument("--rtab", type=Path, default=None,
                    help="--build-store from this finished Rtab instead of calling unitig-caller")
    ap.add_argument("--threads", type=int,
                    default=unitig_cfg.get("threads") or config["preprocessing"]["threads"])
    ap.add_argument("--min-support", type=int, default=int(unitig_cfg["min_support"]))
    args = ap.parse_args()

    print("=" * 78)
    if args.build_store:
        print(f"UNITIG STORE — {args.organism}")
        print("=" * 78)
        build_store(args.organism, config, threads=args.threads,
                    min_support=args.min_support, rtab=args.rtab)
    else:
        print(f"MODEL MATRIX — {args.organism} / {args.antibiotic}")
        print("=" * 78)
        build_model(args.organism, args.antibiotic, config, min_support=args.min_support)


if __name__ == "__main__":
    main()

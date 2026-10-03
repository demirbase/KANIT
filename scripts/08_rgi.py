#!/usr/bin/env python3
"""Step 08 — RGI on every genome of an organism.

Subcommands:

  load     load the CARD release (config card.card_json) into the RGI database
           directory (paths_organism.rgi_db_dir) with `rgi load --local`
  run      `rgi main` on the organism's genomes (contigs, DIAMOND, Perfect and
           Strict hits; Loose hits and nudging stay off), one log per genome; a
           finished genome leaves <genome>.done and is not run again.
           --shard K/N runs every N-th genome starting at K, for a workflow
  collect  every genome's hits into rgi_hits.csv (ORF and protein sequences kept
           for variant, overexpression and rRNA hits only) and rgi_summary.json
           with the RGI and CARD versions; a missing genome stops it

The genomes are those of the organism's unitig store (its panel genomes).
Outputs: paths_organism.rgi_dir.
"""
import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import card_layer  # noqa: E402
from lib.config import get_target, load_config, resolve_path, resolve_tool  # noqa: E402
from lib.databases import card_version  # noqa: E402
from lib.io_utils import run_logged  # noqa: E402
from lib.matrix_store import Store  # noqa: E402

KEEP_SEQUENCES = {"variant", "overexpression", "rrna_variant"}


def _rgi():
    tool = resolve_tool("rgi")
    if not tool:
        sys.exit("ERROR: rgi not found on PATH (or set AMR_RGI_BIN).")
    return tool


def _version(args, cwd):
    out = subprocess.run(args, capture_output=True, text=True, cwd=cwd)
    lines = (out.stdout + out.stderr).strip().splitlines()
    return lines[-1].strip() if lines else None


def genomes_of(organism, config):
    return Store(resolve_path("unitig_store_dir", organism=organism, config=config)).genomes


def load(config):
    db = resolve_path("rgi_db_dir", config=config)
    db.mkdir(parents=True, exist_ok=True)
    card_json = PROJECT_ROOT / config["card"]["card_json"]
    run_logged([_rgi(), "load", "--card_json", card_json, "--local"], db / "rgi_load.log", cwd=db)
    print(f"  ✓ CARD {_version([_rgi(), 'database', '--version', '--local'], db)} loaded in {db}")


def run(organism, config, *, threads, shard):
    db = resolve_path("rgi_db_dir", config=config)
    out = resolve_path("rgi_dir", organism=organism, config=config)
    out.mkdir(parents=True, exist_ok=True)
    fna_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    k, n = shard
    todo = [g for i, g in enumerate(genomes_of(organism, config)) if i % n == k]
    tool = _rgi()
    for i, g in enumerate(todo, 1):
        if (out / f"{g}.done").exists():
            continue
        fna = fna_dir / f"{g}.fna"
        if not fna.exists():
            sys.exit(f"ERROR: assembly missing: {fna}")
        run_logged([tool, "main", "--input_sequence", fna, "--output_file", out / g,
                    "--input_type", "contig", "--alignment_tool", "DIAMOND", "--local",
                    "--clean", "--threads", threads], out / f"{g}.log", cwd=db)
        (out / f"{g}.done").write_text(datetime.now(timezone.utc).isoformat(timespec="seconds"))
        print(f"  ✓ {g} ({i}/{len(todo)})")


def collect(organism, config):
    db = resolve_path("rgi_db_dir", config=config)
    out = resolve_path("rgi_dir", organism=organism, config=config)
    genomes = genomes_of(organism, config)
    missing = [g for g in genomes if not (out / f"{g}.done").exists()]
    if missing:
        sys.exit(f"ERROR: RGI not finished for {len(missing)} genome(s), e.g. {missing[:5]}")
    parts = []
    for g in genomes:
        h = card_layer.read_rgi(out / f"{g}.txt", g)
        drop = ~h["model_type"].isin(KEEP_SEQUENCES)
        h.loc[drop, ["dna", "protein", "ref_protein"]] = ""
        parts.append(h)
    hits = pd.concat(parts, ignore_index=True)
    hits.to_csv(out / "rgi_hits.csv", index=False)
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_genomes": len(genomes), "n_hits": int(len(hits)),
        "hits_by_model_type": hits["model_type"].value_counts().to_dict(),
        "hits_by_cut_off": hits["cut_off"].value_counts().to_dict(),
        "rgi_version": _version([_rgi(), "main", "--version"], db),
        "card_version": _version([_rgi(), "database", "--version", "--local"], db),
        "card_json_version": card_version(PROJECT_ROOT / config["card"]["card_json"]),
    }
    if summary["card_version"] != summary["card_json_version"]:
        sys.exit(f"ERROR: RGI's database is CARD {summary['card_version']}, card.json is "
                 f"{summary['card_json_version']}; run `08_rgi.py load` again.")
    (out / "rgi_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"  ✓ {len(hits):,} hits in {len(genomes):,} genomes -> {out / 'rgi_hits.csv'}")


def _shard(value):
    k, n = (int(x) for x in value.split("/"))
    if not 0 <= k < n:
        raise argparse.ArgumentTypeError("--shard K/N needs 0 <= K < N")
    return k, n


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="RGI on every genome of an organism.")
    ap.add_argument("command", choices=["load", "run", "collect"])
    ap.add_argument("--organism", default=get_target(config=config)[0])
    ap.add_argument("--threads", type=int, default=int(config["preprocessing"]["threads"]))
    ap.add_argument("--shard", type=_shard, default=(0, 1))
    args = ap.parse_args()
    print(f"RGI — {args.organism} — {args.command}")
    if args.command == "load":
        load(config)
    elif args.command == "run":
        run(args.organism, config, threads=args.threads, shard=args.shard)
    else:
        collect(args.organism, config)


if __name__ == "__main__":
    main()

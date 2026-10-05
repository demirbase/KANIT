#!/usr/bin/env python3
"""Step 19 — external validation on CABBAGE (protocol §14, item 4; secondary analysis).

Subcommands:

  download  the release's phenotype table (config cabbage.url) into cabbage_data_dir,
            checked against its MD5; needs internet, so it runs on a login node
  select    the eligible external isolates of every organism (lib/cabbage.py): the rules
            of §2.2 on CABBAGE's fields, without the isolates BV-BRC holds (the genomes.csv
            of every organism's snapshot). Writes to cabbage_dir: cabbage_isolates.csv,
            cabbage_phenotypes.csv, cabbage_pairs.csv and cabbage_selection.json (the
            isolates left after every step)
"""
import argparse
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import cabbage, registry  # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402

USER_AGENT = "KANIT (https://github.com/iumobg/KANIT)"


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def table_path(config: dict) -> Path:
    c = config["cabbage"]
    return resolve_path("cabbage_data_dir", config=config) / str(c["release"]) / "phenotype.parquet"


def download(config: dict) -> Path:
    c = config["cabbage"]
    out = table_path(config)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    req = urllib.request.Request(c["url"], headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
        for block in iter(lambda: r.read(1 << 20), b""):
            f.write(block)
    got = md5_file(tmp)
    if got != c["md5"]:
        tmp.unlink()
        sys.exit(f"ERROR: {c['url']} has MD5 {got}, the protocol names {c['md5']}")
    tmp.replace(out)
    (out.parent / "manifest.json").write_text(json.dumps({
        "release": str(c["release"]), "url": c["url"], "md5": got,
        "bytes": out.stat().st_size,
        "downloaded_at": datetime.now(timezone.utc).isoformat(timespec="seconds")},
        indent=2) + "\n")
    return out


def snapshot_genomes(config: dict) -> pd.DataFrame:
    """genomes.csv of every organism's snapshot: the genomes of its BV-BRC query."""
    parts = []
    for org in registry.load_organisms():
        f = resolve_path("metadata_file", organism=org, config=config).parent / "genomes.csv"
        if not f.exists():
            sys.exit(f"ERROR: no snapshot of {org} ({f}); independence from BV-BRC needs "
                     "every organism's snapshot (run the DOWNLOAD entry)")
        parts.append(pd.read_csv(f, dtype=str, keep_default_na=False,
                                 usecols=["genome_id", "biosample_accession", "sra_accession",
                                          "assembly_accession"]))
    return pd.concat(parts, ignore_index=True)


def select(config: dict) -> dict:
    c = config["cabbage"]
    src = table_path(config)
    if not src.exists():
        sys.exit(f"ERROR: no CABBAGE table at {src}; run `19_cabbage.py download` first")
    md5 = md5_file(src)
    if md5 != c["md5"]:
        sys.exit(f"ERROR: {src} has MD5 {md5}, the protocol names {c['md5']}")
    ph = pd.read_parquet(src, columns=cabbage.COLUMNS)
    ids = cabbage.bvbrc_ids(snapshot_genomes(config))
    isolates, phenotypes, counts = cabbage.select(ph, registry.load_organisms(), ids)
    out = resolve_path("cabbage_dir", config=config)
    out.mkdir(parents=True, exist_ok=True)
    isolates.to_csv(out / "cabbage_isolates.csv", index=False)
    phenotypes.to_csv(out / "cabbage_phenotypes.csv", index=False)
    cabbage.pair_counts(phenotypes).to_csv(out / "cabbage_pairs.csv", index=False)
    summary = {
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "release": str(c["release"]), "md5": md5, "n_records": len(ph),
        "bvbrc_accessions": {k: len(v) for k, v in ids.items()},
        "steps": cabbage.STEPS, "isolates_left": counts,
        "n_isolates": len(isolates), "n_phenotypes": len(phenotypes)}
    (out / "cabbage_selection.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="External validation on CABBAGE.")
    ap.add_argument("command", choices=["download", "select"])
    args = ap.parse_args()
    print(f"CABBAGE — {args.command}")
    if args.command == "download":
        print(f"  ✓ {download(config)}")
    else:
        s = select(config)
        print(f"  {s['n_isolates']} isolates, {s['n_phenotypes']} phenotypes")
        for org, left in s["isolates_left"].items():
            print(f"  {org:24s} " + " → ".join(str(left[k]) for k in cabbage.STEPS))


if __name__ == "__main__":
    main()

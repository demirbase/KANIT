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
  fetch     the GenBank assembly of every eligible isolate from NCBI, whose NCBI taxon is
            at or below the organism's taxa, checked against its genome report (number of
            sequences, length); into cabbage_genomes_dir as <BioSample>.fna. A passed
            assembly already on disk is not fetched again. Writes
            cabbage_download_report.csv and cabbage_fetch.json; needs internet
"""
import argparse
import hashlib
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import bvbrc, cabbage, ncbi_datasets, registry  # noqa: E402
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


REPORT_COLUMNS = ["organism", "biosample_id", "assembly_id", "tax_id", "organism_name", "status",
                  "problem", "sequences", "length", "bytes", "sha256"]


def fetch_one(api, base: dict, contigs: str, length: str, genomes_dir: Path) -> dict:
    try:
        data = api.fasta(base["assembly_id"])
    except Exception as e:                          # recorded; the next fetch tries again
        return {**base, "status": "failed", "problem": f"{type(e).__name__}: {e}"}
    n, total, problem = bvbrc.check_fasta(data, contigs, length)
    if problem:
        return {**base, "status": "failed", "problem": problem, "sequences": n, "length": total}
    genomes_dir.mkdir(parents=True, exist_ok=True)
    path = genomes_dir / f"{base['biosample_id']}.fna"
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)
    return {**base, "status": "passed", "problem": "", "sequences": n, "length": total,
            "bytes": len(data), "sha256": bvbrc.sha256_bytes(data)}


def fetch(config: dict, *, workers: int, api=None) -> pd.DataFrame:
    out = resolve_path("cabbage_dir", config=config)
    iso_file = out / "cabbage_isolates.csv"
    if not iso_file.exists():
        sys.exit(f"ERROR: no {iso_file}; run `19_cabbage.py select` first")
    iso = pd.read_csv(iso_file, dtype=str, keep_default_na=False)
    api = api or ncbi_datasets.Datasets()
    reports = api.genome_reports(iso["assembly_id"])
    numbers = {a: ncbi_datasets.report_numbers(r) for a, r in reports.items()}
    lineages = api.lineages({n[2] for n in numbers.values() if n[2] is not None})
    report_file = out / "cabbage_download_report.csv"
    before = (pd.read_csv(report_file, dtype=str, keep_default_na=False).set_index("biosample_id")
              if report_file.exists() else pd.DataFrame())
    organisms = registry.load_organisms()
    rows, todo = [], []
    for r in iso.itertuples(index=False):
        base = {"organism": r.organism, "biosample_id": r.biosample_id, "assembly_id": r.assembly_id}
        if r.assembly_id not in numbers:
            rows.append({**base, "status": "failed", "problem": "no NCBI genome report"})
            continue
        contigs, length, tax, name = numbers[r.assembly_id]
        base.update(tax_id=tax, organism_name=name)
        if not cabbage.taxon_ok(tax, lineages, organisms[r.organism]["taxids"]):
            rows.append({**base, "status": "excluded",
                         "problem": f"NCBI taxon {tax} ({name}) is not under the organism's taxa"})
            continue
        genomes_dir = resolve_path("cabbage_genomes_dir", organism=r.organism, config=config)
        f = genomes_dir / f"{r.biosample_id}.fna"
        prev = before.loc[r.biosample_id] if r.biosample_id in before.index else None
        if (prev is not None and prev["status"] == "passed" and prev["assembly_id"] == r.assembly_id
                and f.exists() and prev["bytes"] and int(float(prev["bytes"])) == f.stat().st_size):
            rows.append({**prev.to_dict(), "biosample_id": r.biosample_id})
            continue
        todo.append((base, contigs, length, genomes_dir))
    print(f"  {len(iso)} isolates: {len(todo)} to fetch")
    with ThreadPoolExecutor(max_workers=workers) as ex:
        rows += list(ex.map(lambda t: fetch_one(api, *t), todo))
    report = pd.DataFrame(rows).reindex(columns=REPORT_COLUMNS).sort_values(
        ["organism", "biosample_id"], ignore_index=True)
    for c in ("tax_id", "sequences", "length", "bytes"):
        report[c] = pd.to_numeric(report[c], errors="coerce").astype("Int64")
    report.to_csv(report_file, index=False)
    status = report.groupby(["organism", "status"]).size().unstack(fill_value=0)
    (out / "cabbage_fetch.json").write_text(json.dumps({
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_isolates": len(report),
        "status": {o: {k: int(v) for k, v in row.items()} for o, row in status.iterrows()}},
        indent=2) + "\n")
    return report


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="External validation on CABBAGE.")
    ap.add_argument("command", choices=["download", "select", "fetch"])
    ap.add_argument("--workers", type=int, default=4, help="assemblies fetched at once")
    args = ap.parse_args()
    print(f"CABBAGE — {args.command}")
    if args.command == "download":
        print(f"  ✓ {download(config)}")
    elif args.command == "fetch":
        r = fetch(config, workers=args.workers)
        print("  " + ", ".join(f"{k} {v}" for k, v in r["status"].value_counts().items()))
    else:
        s = select(config)
        print(f"  {s['n_isolates']} isolates, {s['n_phenotypes']} phenotypes")
        for org, left in s["isolates_left"].items():
            print(f"  {org:24s} " + " → ".join(str(left[k]) for k in cabbage.STEPS))


if __name__ == "__main__":
    main()

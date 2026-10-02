#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stamp a real Zenodo DOI into the KB (and report what else must change).

Run this ONLY after the deposit exists. It refuses anything that is not a
well-formed Zenodo DOI, because the point of `kb_metadata.zenodo_doi` is that a
copy of the database can identify its own archive — a placeholder in that column
would make verify_artefacts pass while the property it checks is still false,
which is worse than an empty field.

    python scripts/stamp_zenodo_doi.py --db results/kb/kanit.db --doi 10.5281/zenodo.1234567
"""
import argparse
import re
import sqlite3
import sys
from pathlib import Path

DOI_RE = re.compile(r"^10\.5281/zenodo\.\d{4,}$")
# Every place the thesis names the DOI. Kept here so the stamp and the prose cannot
# drift: the script does not edit the thesis (it is not in this repo), it lists them.
THESIS_SPOTS = [
    "3_malzeme_yontem_EN.md §3.8.2 — replace 10.5281/zenodo.XXXXXXX with the concept DOI",
    "3_malzeme_yontem_EN.md §3.8.3 — F1 and F3 move from 'not met' to met; rewrite both entries",
    "3_malzeme_yontem_EN.md §3.8.2 — drop 'In the delivered database that column is empty'",
    "4.1_kb_genel_gorunum_EN.md — the SAYI KAYNAKLARI note about the NULL column",
    "5.6_acik_bilim_surdurulebilirlik_EN.md — the open-science section's DOI statement",
    "scripts/kb_fair_mapping.py — F1/F3 status, then regenerate fair_mapping.csv",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="results/kb/kanit.db")
    ap.add_argument("--doi", required=True, help="concept DOI, e.g. 10.5281/zenodo.1234567")
    a = ap.parse_args()

    doi = a.doi.strip().replace("https://doi.org/", "")
    if not DOI_RE.match(doi):
        sys.exit(f"ERROR: {doi!r} is not a Zenodo DOI of the form 10.5281/zenodo.<digits>. "
                 "Refusing to stamp a placeholder — an empty column is honest, a fake one is not.")

    conn = sqlite3.connect(a.db)
    cur = conn.execute("SELECT zenodo_doi FROM kb_metadata").fetchone()
    print(f"KB: {a.db}\n  zenodo_doi now: {cur[0] or '(empty)'}\n  stamping: {doi}")
    with conn:
        conn.execute("UPDATE kb_metadata SET zenodo_doi = ?", (doi,))
    got = conn.execute("SELECT zenodo_doi FROM kb_metadata").fetchone()[0]
    conn.close()
    print(f"  ✓ stamped: {got}\n")
    print("Still to change by hand — the DOI appears in prose too:")
    for s in THESIS_SPOTS:
        print(f"  - {s}")
    print("\nThen: python scripts/verify_artefacts.py   (expects 45/45)")


if __name__ == "__main__":
    main()

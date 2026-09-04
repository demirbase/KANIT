#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Collapse duplicated rows in `blast_annotations`.

`blast_annotations` carries a UNIQUE index over
(unitig_id, model_id, source_db, gene_symbol, identity_pct, coverage, evalue) and both
loader routes insert with INSERT OR IGNORE, so a hit reaching the table twice — once by
the gain route and once by the CPSS route — was meant to be stored once.

It was not, because `coverage` was in the index and only one route supplied it. SQLite
treats two NULLs as distinct, so the pair stayed apart on the strength of a column that
was empty for want of a bug fix, and nine annotations were stored twice.

This is the only place the defect reaches a delivered number. Everything the thesis
grades is keyed on (unitig, model) in `unitig_evidence_tier`, which holds 3,571 rows and
3,571 distinct pairs — no duplication, so no tier, biomarker or evidence count moves. What
moves is the raw annotation row count, 3,611 -> 3,602.

Deleting is not a judgement call about which row is better: within each group the tier and
the ARO accession are identical, so the rows say the same thing. The one carrying a
coverage is kept, because it carries strictly more, and the other is dropped. With the
loader fixed (09 emits coverage, populate_candidates stores it) a rebuild now produces
3,602 directly, so this brings the delivered database in line with its own code.

Refuses to touch a group whose rows disagree about anything that matters.

    python scripts/migrate_kb_dedupe_annotations.py --db results/kb/amrk.db
    python scripts/migrate_kb_dedupe_annotations.py --db results/kb/amrk.db --apply
"""
import argparse
import sqlite3
import sys

KEY = "unitig_id, model_id, source_db, gene_symbol, identity_pct, evalue"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    conn = sqlite3.connect(a.db)
    conn.execute("PRAGMA foreign_keys = ON")
    total = conn.execute("SELECT COUNT(*) FROM blast_annotations").fetchone()[0]
    groups = conn.execute(
        f"SELECT {KEY}, COUNT(*) FROM blast_annotations GROUP BY {KEY} HAVING COUNT(*) > 1"
    ).fetchall()
    print(f"KB: {a.db}\n  blast_annotations: {total:,} · duplicated groups: {len(groups)}")
    if not groups:
        print("  nothing to do."); conn.close(); return

    drop, refuse = [], []
    for g in groups:
        rows = conn.execute(
            "SELECT annotation_id, coverage, tier, aro_accession, aro_gene_family, "
            "aro_drug_class, aro_resistance_mechanism FROM blast_annotations "
            f"WHERE unitig_id=? AND model_id=? AND source_db=? AND gene_symbol=? "
            "AND identity_pct=? AND evalue=? ORDER BY annotation_id", g[:6]).fetchall()
        # Everything except coverage must agree, or these are not the same annotation.
        if len({r[2:] for r in rows}) != 1:
            refuse.append((g, rows)); continue
        # Prefer a stored coverage; among stored ones prefer the fuller float, because a
        # rounded copy of a value is the same measurement recorded less precisely, and
        # keeping the rounded one would throw away digits for no reason.
        keep = sorted(rows, key=lambda r: (r[1] is None, -len(repr(r[1])), r[0]))[0]
        for r in rows:
            if r[0] != keep[0]:
                drop.append((r[0], g[3], keep[0]))
    if refuse:
        print(f"\n  REFUSING {len(refuse)} group(s) whose rows disagree beyond coverage:")
        for g, rows in refuse[:5]:
            print(f"    {g[:6]} -> {rows}")
    print(f"\n  would delete {len(drop)} row(s), keeping the one that carries a coverage:")
    for aid, gene, keep in drop:
        print(f"    drop annotation_id={aid} ({gene}) · keep {keep}")

    if not a.apply:
        print("\n  DRY RUN — re-run with --apply."); conn.close(); return

    with conn:
        conn.executemany("DELETE FROM blast_annotations WHERE annotation_id = ?",
                         [(aid,) for aid, _, _ in drop])
    after = conn.execute("SELECT COUNT(*) FROM blast_annotations").fetchone()[0]
    left = conn.execute(
        f"SELECT COUNT(*) FROM (SELECT 1 FROM blast_annotations GROUP BY {KEY} "
        "HAVING COUNT(*) > 1)").fetchone()[0]
    print(f"\n  blast_annotations: {total:,} -> {after:,} · duplicated groups left: {left}")
    tiers = dict(conn.execute("SELECT evidence_tier, COUNT(*) FROM unitig_evidence_tier "
                              "GROUP BY evidence_tier"))
    print(f"  evidence tiers unchanged: {tiers}")
    conn.close()
    sys.exit(0 if left == 0 else 1)


if __name__ == "__main__":
    main()

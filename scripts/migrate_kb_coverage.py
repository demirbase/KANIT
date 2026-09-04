#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Backfill `blast_annotations.coverage` for the rows one loader path never wrote.

populate_database inserts BLAST annotations by two routes. The CPSS route
(populate_cpss) lists `coverage` in its column list and 13b supplies the value.
The gain route (populate_candidates) omits the column from its INSERT entirely, so
every annotation arriving that way stores NULL — 416 rows carrying a real tier.

Grading was never affected: step 09 computes coverage from the alignment and applies
it to the tier cutoffs before anything is written. Only the stored number is absent.
So this recomputes it the way 09 does, `length / qlen`, from the delivered CARD BLAST
TSVs, and it does not trust the recomputation blindly:

  * the source row is found by an EXACT match on the (identity, E-value) pair the KB
    already stores, not by re-deriving 09's best-hit rule — if the KB and the TSV
    disagree the row is skipped rather than guessed;
  * a value above 1.0 is impossible for an alignment against its own query and is
    treated as a mismatch;
  * every recovered value is checked against the floor its recorded tier requires
    (confirmed 0.95, candidate 0.80, weak 0.60) and the migration aborts if any row
    fails, because that would mean the row and the hit are not the same thing.

Rows whose gene_symbol is the literal string "nan" are left alone: they hold no
identity, E-value, ARO or hit of any kind and are not annotations (see the note in
the repository's audit log). Two real rows could not be matched and stay NULL.

Safe to re-run: rows that already carry a coverage are never touched.

    python scripts/migrate_kb_coverage.py --db results/kb/amrk.db
    python scripts/migrate_kb_coverage.py --db results/kb/amrk.db --apply
"""
import argparse
import glob
import os
import sqlite3
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent

COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
        'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'stitle']
FLOOR = {"confirmed": 0.95, "candidate": 0.80, "weak": 0.60, "none": 0.0}


def recover(conn, results_root):
    kb = pd.read_sql("""
        SELECT ba.annotation_id, ba.tier, u.sequence, ba.identity_pct, ba.evalue,
               pr.organism, m.antibiotic
        FROM blast_annotations ba
        JOIN unitigs u        ON u.unitig_id = ba.unitig_id
        JOIN models m         ON m.model_id  = ba.model_id
        JOIN pipeline_runs pr ON pr.run_id   = m.run_id
        WHERE ba.coverage IS NULL AND ba.gene_symbol <> 'nan'""", conn)
    out, skipped = {}, {"no files": 0, "not in FASTA": 0, "no matching hit": 0,
                        "qlen missing": 0, "coverage > 1": 0}
    for (org, ab), grp in kb.groupby(["organism", "antibiotic"]):
        d = Path(results_root) / org / ab / "05_explainability"
        fa = glob.glob(str(d / f"02_top_*_features_{ab}.fasta"))
        ts = d / f"03_card_blast_results_{ab}.tsv"
        if not fa or not ts.exists():
            skipped["no files"] += len(grp); continue
        seq2id, hdr = {}, None
        with open(fa[0], encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    hdr = line[1:]
                elif hdr:
                    seq2id[line] = hdr; hdr = None
        df = pd.read_csv(ts, sep="\t", header=None, names=COLS, dtype=str)
        for c in ("pident", "evalue", "length", "qlen"):
            df[c] = pd.to_numeric(df[c], errors="coerce")
        for r in grp.itertuples():
            q = seq2id.get(r.sequence)
            if q is None:
                skipped["not in FASTA"] += 1; continue
            if pd.isna(r.identity_pct) or pd.isna(r.evalue):
                skipped["no matching hit"] += 1; continue
            h = df[(df.qseqid == q)
                   & (df.pident.round(3) == round(float(r.identity_pct), 3))
                   & (df.evalue == float(r.evalue))]
            if h.empty:
                skipped["no matching hit"] += 1; continue
            b = h.iloc[0]
            if not (b.qlen > 0):
                skipped["qlen missing"] += 1; continue
            cov = float(b.length) / float(b.qlen)
            if cov > 1.0:
                skipped["coverage > 1"] += 1; continue
            out[int(r.annotation_id)] = (round(cov, 6), r.tier)
    return len(kb), out, skipped


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--results", default="results")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    conn = sqlite3.connect(a.db)
    conn.execute("PRAGMA foreign_keys = ON")
    before = conn.execute(
        "SELECT COUNT(*), SUM(coverage IS NOT NULL) FROM blast_annotations").fetchone()
    print(f"KB: {a.db}\n  blast_annotations: {before[0]:,} · coverage set: {before[1]:,}")

    target, rec, skipped = recover(conn, a.results)
    print(f"  real annotations missing coverage: {target} · recovered: {len(rec)}")
    for k, v in skipped.items():
        if v:
            print(f"    skipped, {k}: {v}")

    bad = [(aid, t, c) for aid, (c, t) in rec.items() if c < FLOOR.get(t, 0.0)]
    if bad:
        conn.close()
        sys.exit(f"ERROR: {len(bad)} recovered value(s) fall below the floor their recorded "
                 f"tier requires, e.g. {bad[:3]}. The row and the hit are not the same "
                 "alignment — refusing to write.")
    print(f"  ✓ all {len(rec)} recovered values satisfy their tier's coverage floor")

    if not rec:
        conn.close(); return
    # Filling the value can collide with the table's UNIQUE index, because SQLite treats
    # two NULLs as distinct: where the same hit was loaded by BOTH routes, the rows stayed
    # apart only for want of a coverage. Supplying it makes them identical. Those are real
    # duplicates, but removing them would change a row count the thesis and the artefact
    # check both pin, so they are reported and left NULL — a separate decision, not a
    # side effect of this one.
    dup = []
    for aid, (cov, tier) in list(rec.items()):
        r = conn.execute("""SELECT unitig_id, model_id, source_db, gene_symbol,
                                   identity_pct, evalue FROM blast_annotations
                            WHERE annotation_id = ?""", (aid,)).fetchone()
        clash = conn.execute("""SELECT COUNT(*) FROM blast_annotations
              WHERE unitig_id=? AND model_id=? AND source_db=? AND gene_symbol=?
                AND identity_pct=? AND evalue=? AND coverage=? AND annotation_id<>?""",
                             (*r, cov, aid)).fetchone()[0]
        if clash:
            dup.append(aid); rec.pop(aid)
    if dup:
        print(f"  {len(dup)} row(s) left NULL: filling them would duplicate a row that "
              "already carries the same hit and coverage (loaded by both routes)")

    if not a.apply:
        print("\n  DRY RUN — re-run with --apply.")
        conn.close(); return

    with conn:
        conn.executemany("UPDATE blast_annotations SET coverage = ? WHERE annotation_id = ?",
                         [(c, aid) for aid, (c, t) in rec.items()])
    after = conn.execute(
        "SELECT COUNT(*), SUM(coverage IS NOT NULL) FROM blast_annotations").fetchone()
    real = conn.execute(
        "SELECT COUNT(*), SUM(coverage IS NOT NULL) FROM blast_annotations "
        "WHERE gene_symbol <> 'nan'").fetchone()
    print(f"\n  coverage set: {before[1]:,} -> {after[1]:,}")
    print(f"  of the {real[0]:,} real annotations, {real[1]:,} now carry a coverage "
          f"({100*real[1]/real[0]:.1f}%)")
    conn.close()


if __name__ == "__main__":
    main()

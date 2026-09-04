#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Correct `unitigs.k` on a populated AMRK-DB: the k-mer baseline's k was stamped
on unitig rows.

populate_database read `preprocessing.k_length` unconditionally when inserting into
`unitigs`. That key is the KMC baseline's k (21). The delivered run uses the unitig
representation, whose de Bruijn graph is built at 31 — 03u passes no `--kmer`, so
unitig-caller's own default stands. Every unitig row therefore records a k the graph
was not built at.

Confirmed three independent ways before this migration was written:
  * unitig_caller/__main__.py in the delivered container: `--kmer` default = 31
  * 03u_unitig_matrix.py invokes `--call --refs --rtab --out --threads`, no --kmer
  * the delivered store contains no unitig shorter than 31 bp, which a graph at
    k = 31 cannot produce and one at k = 21 could not avoid

Nothing else in the KB derives from this column: BLAST coverage is computed against
each unitig's own qlen, never against k. This is a metadata correction, not a
recomputation, which is why it is a migration rather than a rebuild — rebuilding
would put the only repaired copy of the KB back through a loader for no gain.

The loader itself is fixed in the same change (populate_database._feature_k), so a
future rebuild does not reintroduce this.

Safe to re-run: it is a no-op once every row is at the target k.

    python scripts/migrate_kb_unitig_k.py --db results/kb/amrk.db
    python scripts/migrate_kb_unitig_k.py --db results/kb/amrk.db --apply
"""
import argparse
import sqlite3
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--apply", action="store_true",
                    help="write the change (default: report what would change)")
    a = ap.parse_args()

    target = (load_config().get("unitig", {}) or {}).get("k")
    if target is None:
        sys.exit("ERROR: config has no unitig.k — nothing to migrate to.")
    target = int(target)

    conn = sqlite3.connect(a.db)
    conn.execute("PRAGMA foreign_keys = ON")
    before = dict(conn.execute("SELECT k, COUNT(*) FROM unitigs GROUP BY k"))
    total = sum(before.values())
    print(f"KB: {a.db}")
    print(f"  unitigs: {total:,} · k distribution now: {before}")
    print(f"  target k (config unitig.k): {target}")

    # The shortest sequence is a direct check on the target: no path through a
    # de Bruijn graph can be shorter than one node, so min(length) >= k, and for a
    # graph that has any tip at all it is equal to it.
    shortest = conn.execute("SELECT MIN(LENGTH(sequence)) FROM unitigs").fetchone()[0]
    print(f"  shortest unitig in the KB: {shortest} bp", end="")
    if shortest < target:
        conn.close()
        sys.exit(f"\nERROR: a unitig is shorter than the target k ({shortest} < {target}). "
                 "That is impossible for a graph built at this k — refusing to write.")
    print(" (>= k, consistent)" if shortest > target else " (== k, as expected)")

    stale = int(conn.execute("SELECT COUNT(*) FROM unitigs WHERE k != ?", (target,)).fetchone()[0])
    if stale == 0:
        print("  nothing to do — every row already at the target k.")
        conn.close()
        return

    if not a.apply:
        print(f"\n  DRY RUN — {stale:,} row(s) would move to k = {target}. Re-run with --apply.")
        conn.close()
        return

    with conn:
        conn.execute("UPDATE unitigs SET k = ?", (target,))
    after = dict(conn.execute("SELECT k, COUNT(*) FROM unitigs GROUP BY k"))
    ok = after == {target: total}
    print(f"\n  updated {stale:,} row(s) · k distribution now: {after}")
    print("  ✓ every unitig at the target k" if ok else "  ✗ UNEXPECTED distribution")
    conn.close()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

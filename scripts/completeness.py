#!/usr/bin/env python3
"""Completeness of a workflow run (plan F3.6): every expected output of every step, for
every organism and panel model, exists and is valid.

Valid: a JSON file parses; a CSV file parses, has the required columns and at least
the required rows; any other file is not empty. A step a model does not reach is "n/a":
the steps after the cross-validation design for a model that is not evaluable.

  completeness.py check --entry main --organisms ecoli,kpneumoniae --out-dir runs/nextflow/<run>
  completeness.py gate --out-dir runs/nextflow/<run>

Writes completeness.csv (one row per check), completeness_matrix.csv (unit × step:
ok, missing, invalid or n/a) and completeness_summary.json ("complete": every check
ok or n/a). The run's last task fails when the run is not complete.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import panel  # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402


@dataclass(frozen=True)
class Check:
    step: str
    key: str                 # paths_organism key; a *_file key is the file itself or its directory
    name: str                # file under the key's directory ("" = the *_file itself)
    kind: str = "csv"        # csv, json or file
    columns: tuple = ()
    min_rows: int = 0
    evaluable_only: bool = False


SNAPSHOT = [Check("snapshot", "metadata_file", "snapshot.json", "json"),
            Check("snapshot", "metadata_file", "", "csv", ("Genome ID",), 1),
            Check("snapshot", "metadata_file", "genomes.csv", "csv", ("genome_id",), 1)]
ORGANISM = SNAPSHOT + [
    Check("genome_qc", "genome_qc_dir", "02d_genome_qc_{organism}.csv", "csv",
          ("genome_id", "pass_overall"), 1),
    Check("genome_qc", "genome_qc_dir", "versions.json", "json"),
    Check("lineage", "lineage_dir", "poppunk_clusters.csv", "csv", ("Genome ID", "Cluster"), 1),
    Check("lineage", "lineage_dir", "versions.json", "json")]
MODELLED_ORGANISM = [
    Check("unitig_store", "unitig_store_dir", "store_summary.json", "json"),
    Check("rgi", "rgi_dir", "rgi_hits.csv", "csv", ("genome_id", "aro")),
    Check("rgi", "rgi_dir", "rgi_summary.json", "json"),
    Check("external", "external_dir", "comparison.csv", "csv", ("model_id", "tool"), 1),
    Check("external", "external_dir", "amrfinder_calls.csv", "csv", ("genome_id",)),
    Check("external", "external_dir", "resfinder_calls.csv", "csv", ("genome_id", "antibiotic")),
    Check("external", "external_dir", "versions.json", "json")]
MODEL = [
    Check("matrix", "matrix_dir", "matrix_summary.json", "json"),
    Check("matrix", "matrix_dir", "genomes.csv", "csv", ("Genome ID", "label"), 1),
    Check("matrix", "matrix_dir", "patterns.csv", "csv", ("pattern_id",)),
    Check("cv", "cv_dir", "cv_design.json", "json"),
    Check("cv", "cv_dir", "metrics.json", "json", evaluable_only=True),
    Check("cv", "cv_dir", "oof_predictions.csv", "csv", ("genome_id", "p"), 1, True),
    Check("cv", "cv_dir", "final/record.json", "json", evaluable_only=True),
    Check("label_permutation", "label_permutation_dir", "label_permutation.json", "json",
          evaluable_only=True),
    Check("cpss", "cpss_dir", "cpss.csv", "csv", ("pattern_id", "pi"), evaluable_only=True),
    Check("cpss", "candidates_file", "", "csv", ("pattern_id",), evaluable_only=True),
    Check("prevalence", "layers_dir", "prevalence.csv", "csv", ("pattern_id", "passes"),
          evaluable_only=True),
    Check("mda", "layers_dir", "mda.csv", "csv", ("pattern_id", "passes"), evaluable_only=True),
    Check("mda", "layers_dir", "mda_clusters.csv", "csv", evaluable_only=True),
    Check("pyseer", "pyseer_dir", "pyseer_tested.csv", "csv", evaluable_only=True),
    Check("pyseer", "pyseer_dir", "pyseer_summary.json", "json", evaluable_only=True),
    Check("card", "card_layer_dir", "card_unitigs.csv", "csv", ("pattern_id",),
          evaluable_only=True),
    Check("grading", "grades_dir", "grades_patterns.csv", "csv", ("pattern_id",),
          evaluable_only=True),
    Check("grading", "grades_dir", "grades_summary.json", "json", evaluable_only=True)]
CONTEXT = [Check("context", "context_dir", "unitig_context.csv", "csv", ("unitig_id",))]
GLOBAL_MAIN = [Check("panel", "panel_dir", "panel_decisions.csv", "csv",
                     ("organism", "antibiotic", "decision"), 1)]
ACROSS = [Check("label_permutation", "cross_model_dir", "label_permutation.csv", "csv")]
KB = [Check("kb", "kb_dir", "kanit.sqlite", "file"),
      Check("kb", "kb_dir", "kb_report.json", "json"),
      Check("hypotheses", "cross_model_dir", "hypotheses/hypotheses.json", "json")]


def path_of(c: Check, config: dict, organism=None, antibiotic=None) -> Path:
    p = resolve_path(c.key, organism=organism, antibiotic=antibiotic, config=config)
    if c.key.endswith("_file"):
        return p if not c.name else p.parent / c.name
    return p / c.name.format(organism=organism)


def validate(c: Check, p: Path) -> tuple[str, str]:
    if not p.exists():
        return "missing", ""
    if p.stat().st_size == 0:
        return "invalid", "empty file"
    try:
        if c.kind == "json":
            json.loads(p.read_text())
        elif c.kind == "csv":
            df = pd.read_csv(p, nrows=None if c.min_rows else 0)
            lacking = [x for x in c.columns if x not in df.columns]
            if lacking:
                return "invalid", f"no column {lacking}"
            if len(df) < c.min_rows:
                return "invalid", f"{len(df)} rows, needs {c.min_rows}"
    except (ValueError, pd.errors.ParserError, UnicodeDecodeError) as e:
        return "invalid", f"{type(e).__name__}: {e}"
    return "ok", ""


def run_checks(entry: str, organisms: list[str], config: dict) -> pd.DataFrame:
    rows: list[dict] = []

    def check(unit, checks, *, organism=None, antibiotic=None, evaluable=True):
        for c in checks:
            p = path_of(c, config, organism, antibiotic)
            status, detail = (("n/a", "not evaluable") if c.evaluable_only and not evaluable
                              else validate(c, p))
            rows.append({"unit": unit, "step": c.step, "file": str(p), "status": status,
                         "detail": detail})

    if entry == "DOWNLOAD":
        for org in organisms:
            check(org, SNAPSHOT, organism=org)
        return pd.DataFrame(rows)
    if entry == "KB":
        check("all", KB)
        return pd.DataFrame(rows)
    decisions_file = resolve_path("panel_dir", config=config) / "panel_decisions.csv"
    decisions = (panel.read_panel(decisions_file) if decisions_file.exists()
                 else pd.DataFrame(columns=["organism", "antibiotic", "decision"]))
    included = decisions[(decisions["decision"] == panel.INCLUDED)
                         & decisions["organism"].isin(organisms)]
    modelled = sorted(set(included["organism"]))
    if entry == "CONTEXT":
        for org in modelled:
            check(org, CONTEXT, organism=org)
        return pd.DataFrame(rows)
    check("all", GLOBAL_MAIN)
    for org in organisms:
        check(org, ORGANISM, organism=org)
    for org in modelled:
        check(org, MODELLED_ORGANISM, organism=org)
    any_evaluable = False
    for row in included.itertuples():
        design = resolve_path("cv_dir", organism=row.organism, antibiotic=row.antibiotic,
                              config=config) / "cv_design.json"
        evaluable = (not design.exists()) or bool(json.loads(design.read_text()).get("evaluable"))
        any_evaluable |= evaluable
        check(f"{row.organism}__{row.antibiotic}", MODEL, organism=row.organism,
              antibiotic=row.antibiotic, evaluable=evaluable)
    if any_evaluable:
        check("all", ACROSS)
    return pd.DataFrame(rows)


def matrix(table: pd.DataFrame) -> pd.DataFrame:
    order = {"missing": 3, "invalid": 2, "ok": 1, "n/a": 0}
    worst = (table.assign(rank=table["status"].map(order))
             .groupby(["unit", "step"], sort=False)["rank"].max().map({v: k for k, v in
                                                                      order.items()}))
    steps = list(dict.fromkeys(table["step"]))
    return worst.unstack("step").reindex(columns=steps).fillna("").reset_index()


def gate(out_dir: Path) -> None:
    """Fail when the run is not complete (the workflow's last task)."""
    s = json.loads((out_dir / "completeness_summary.json").read_text())
    if not s["complete"]:
        lines = [f"  {f['unit']} {f['step']}: {f['status']} {f['file']} {f['detail']}"
                 for f in s["failures"][:20]]
        sys.exit(f"INCOMPLETE RUN — {len(s['failures'])} missing or invalid output(s):\n"
                 + "\n".join(lines))
    print(f"COMPLETE — {s['percent_ok']}% of {s['applicable']} checks ok")


def main():
    ap = argparse.ArgumentParser(description="Completeness of a workflow run.")
    ap.add_argument("command", choices=["check", "gate"])
    ap.add_argument("--entry", choices=["main", "DOWNLOAD", "CONTEXT", "KB"])
    ap.add_argument("--organisms", help="comma-separated registry ids")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()
    if args.command == "gate":
        gate(args.out_dir)
        return
    if not args.entry or not args.organisms:
        ap.error("check needs --entry and --organisms")
    table = run_checks(args.entry, args.organisms.split(","), load_config())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out_dir / "completeness.csv", index=False)
    if not table.empty:
        matrix(table).to_csv(args.out_dir / "completeness_matrix.csv", index=False)
    bad = table[table["status"].isin(["missing", "invalid"])]
    applicable = table[table["status"] != "n/a"]
    summary = {"entry": args.entry, "complete": bad.empty, "checks": len(table),
               "applicable": len(applicable),
               "percent_ok": round(100 * (applicable["status"] == "ok").mean(), 2)
               if len(applicable) else 100.0,
               "failures": bad[["unit", "step", "file", "status", "detail"]].to_dict("records")}
    (args.out_dir / "completeness_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"COMPLETENESS — {args.entry}: {summary['percent_ok']}% of {summary['applicable']} "
          f"checks ok; {len(bad)} missing or invalid")


if __name__ == "__main__":
    main()

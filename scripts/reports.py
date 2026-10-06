#!/usr/bin/env python3
"""Figures, model reports and the numbers file (plan F4.2), from the steps' CSV tables.

  reports.py --entry main --organisms ecoli,kpneumoniae

At the end of every entry it draws the QC figure of every step whose tables exist
(lib/figures.py; Appendix C), writes for every panel model a self-contained HTML
report (data, cross-validation, label permutation, evidence layers, grades, the best
candidates, the comparison with genotype-based tools and the model's figures), and
rewrites tez_sayilari.csv: every number the thesis and the article quote, with the
table and column it comes from. A figure whose table is missing is skipped and named
in reports_summary.json (completeness.py reports the missing table).

Outputs: paths_organism.reports_dir (global and organism figures, the numbers file)
and model_reports_dir (a model's figures and report.html).
"""
from __future__ import annotations

import argparse
import base64
import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import contract, panel  # noqa: E402
from lib import figures as F  # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402
from lib.run_metadata import git_commit_hash  # noqa: E402

GRADE_ORDER = {g: i for i, g in enumerate(F.GRADES)}


class Tables:
    """The contract's tables of one run, read when asked for (None when absent)."""

    def __init__(self, config: dict):
        self.config = config
        self.skipped: list[str] = []

    def get(self, tid: str, organism=None, antibiotic=None) -> pd.DataFrame | None:
        try:
            p = contract.table_path(tid, self.config, organism, antibiotic)
        except KeyError:
            return None
        if not p.exists():
            return None
        return pd.read_csv(p, keep_default_na=True, low_memory=False)

    def json(self, location, file, organism=None, antibiotic=None) -> dict | None:
        p = contract.output_path(location, file, self.config, organism, antibiotic)
        return json.loads(p.read_text()) if p.exists() else None


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    F.plt.close(fig)
    return path


# ---- figures -----------------------------------------------------------------------
def organism_figures(t: Tables, org: str, out: Path) -> list[Path]:
    done = []
    long = t.get("phenotypes_long", org)
    if long is not None:
        done.append(_save(F.phenotype_counts(long, f"{org} — phenotypes"), out / "snapshot.png"))
    qc = t.get("genome_qc", org)
    if qc is not None:
        summ = t.json("genome_qc_dir", "02d_genome_qc_summary_{organism}.json", org) or {}
        done.append(_save(F.genome_qc(qc, summ.get("thresholds", {}), f"{org} — genome QC"),
                          out / "genome_qc.png"))
    cl = t.get("lineage_clusters", org)
    if cl is not None:
        done.append(_save(F.lineage_sizes(cl, f"{org} — lineages"), out / "lineages.png"))
    hits = t.get("rgi_hits", org)
    if hits is not None:
        done.append(_save(F.rgi_hits(hits, f"{org} — RGI"), out / "rgi.png"))
    comp = t.get("external_comparison", org)
    if comp is not None:
        done.append(_save(F.external_agreement(comp, f"{org} — genotype-based prediction"),
                          out / "external.png"))
    ctx = t.get("unitig_context", org)
    if ctx is not None:
        done.append(_save(F.context_plasmid(ctx, f"{org} — NCBI context"), out / "context.png"))
    return done


def model_figures(t: Tables, org: str, ab: str, out: Path) -> list[tuple[str, Path]]:
    title = f"{org} / {ab}"
    done: list[tuple[str, Path]] = []

    def add(name, fig):
        done.append((name, _save(fig, out / f"{name}.png")))

    pats, mg = t.get("model_patterns", org, ab), t.get("model_genomes", org, ab)
    if pats is not None and mg is not None:
        add("matrix", F.matrix_patterns(pats, len(mg), title))
    oof = t.get("oof_predictions", org, ab)
    rel, fm = t.get("reliability", org, ab), t.get("fold_metrics", org, ab)
    comp = t.get("fold_composition", org, ab)
    if all(x is not None for x in (oof, rel, fm, comp)):
        add("cross_validation", F.cv_quality(oof, rel, fm, comp, title))
    null, lp = t.get("label_permutation_null", org, ab), t.get("label_permutation_models")
    if null is not None:
        row = (lp[(lp["organism"] == org) & (lp["antibiotic"] == ab)]
               if lp is not None else pd.DataFrame())
        obs = float(row["auc_observed"].iloc[0]) if len(row) else None
        add("label_permutation", F.permutation_null(null, obs, title))
    cpss = t.get("cpss_patterns", org, ab)
    if cpss is not None:
        add("cpss", F.cpss_frequency(cpss, float(t.config["cpss"]["pi_threshold"]), title))
    prev = t.get("prevalence_layer", org, ab)
    if prev is not None:
        add("prevalence", F.prevalence_volcano(prev, title))
    mda = t.get("mda_layer", org, ab)
    if mda is not None:
        add("mda", F.mda_scatter(mda, title))
    qq = t.get("pyseer_qq", org, ab)
    if qq is not None:
        s = t.json("pyseer_dir", "pyseer_summary.json", org, ab) or {}
        lam = {"tested": s.get("lambda_tested"), "background": s.get("lambda_background")}
        add("pyseer", F.pyseer_qq(qq, title, lam))
    card = t.get("card_patterns", org, ab)
    if card is not None:
        add("card", F.card_states(card, title))
    grades = t.get("grades_patterns", org, ab)
    if grades is not None:
        add("grading", F.evidence_funnel(grades, title))
    return done


def global_figures(t: Tables, out: Path) -> list[Path]:
    done = []
    dec = t.get("panel_decisions")
    if dec is not None:
        done.append(_save(F.panel_minority(dec, int(t.config["panel"]["min_minority"]), "Panel"),
                          out / "panel.png"))
    lp = t.get("label_permutation_models")
    if lp is not None:
        done.append(_save(F.permutation_models(lp, "Label permutation"),
                          out / "label_permutation.png"))
    h3 = t.get("h3_pairs")
    if h3 is not None:
        for org in sorted(h3["organism_id"].unique()):
            done.append(_save(F.overlap_heatmap(h3, org), out / org / "candidate_overlap.png"))
    return done


# ---- model report ------------------------------------------------------------------
CSS = """body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;
color:#222}h1{margin-bottom:0}.meta{color:#666}table{border-collapse:collapse;margin:.5em 0 1.5em;
font-size:13px}th,td{border:1px solid #ddd;padding:3px 8px;text-align:right}th{background:#f4f4f4}
td:first-child,th:first-child{text-align:left}img{max-width:100%;margin:.5em 0 1.5em}"""


def _table(df: pd.DataFrame | None) -> str:
    if df is None or df.empty:
        return "<p class='meta'>not available</p>"
    return df.to_html(index=False, float_format=lambda x: f"{x:.4g}", border=0, na_rep="")


def model_report(t: Tables, org: str, ab: str, figs: list[tuple[str, Path]], out: Path) -> Path:
    mg = t.get("model_genomes", org, ab)
    data = None
    if mg is not None:
        data = pd.DataFrame([{"genomes": len(mg), "resistant": int((mg["label"] == 1).sum()),
                              "susceptible": int((mg["label"] == 0).sum()),
                              "lineages": mg["lineage"].nunique()}])
    rm = t.get("repeat_metrics", org, ab)
    cv = None
    if rm is not None:
        cols = [c for c in ("roc_auc", "pr_auc", "balanced_accuracy", "mcc", "brier") if c in rm]
        g = rm.groupby("arm")[cols]
        cv = (g.mean().round(4).astype(str) + " ± " + g.std().round(4).astype(str)).reset_index()
    lp = t.get("label_permutation_models")
    lp_row = lp[(lp["organism"] == org) & (lp["antibiotic"] == ab)] if lp is not None else None
    grades = t.get("grades_patterns", org, ab)
    layers = counts = best = None
    if grades is not None:
        layers = pd.DataFrame([{k: int(F._bool(grades[k]).sum())
                                for k in ("prevalence", "mda", "cpss", "pyseer")}])
        counts = pd.DataFrame({m: grades[f"{m}_grade"].value_counts().reindex(F.GRADES,
                                                                            fill_value=0)
                               for m in ("allele_aware", "homolog_only")}).reset_index(
            names="grade")
        best = (grades.assign(_o=grades["grade"].map(GRADE_ORDER))
                .sort_values(["_o", "n_layers"], ascending=[True, False]).head(25)
                [["pattern_id", "n_members", "layers_passed", "allele_aware_card_state",
                  "allele_aware_grade", "homolog_only_grade"]])
    comp = t.get("external_comparison", org)
    ext = comp[comp["model_id"] == f"{org}__{ab}"] if comp is not None else None
    parts = [f"<h1>{html.escape(org)} / {html.escape(ab)}</h1>",
             f"<p class='meta'>KANIT v1.0 · {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC · "
             f"commit {html.escape(str(git_commit_hash() or 'unknown'))}</p>",
             "<h2>Data</h2>", _table(data),
             "<h2>Cross-validation (mean ± SD over repeats)</h2>", _table(cv),
             "<h2>Label permutation</h2>", _table(lp_row),
             "<h2>Evidence layers that fire</h2>", _table(layers),
             "<h2>Grades</h2>", _table(counts),
             "<h2>Best candidates</h2>", _table(best),
             "<h2>Genotype-based prediction</h2>", _table(ext), "<h2>Figures</h2>"]
    for name, path in figs:
        b64 = base64.b64encode(path.read_bytes()).decode()
        parts.append(f"<h3>{html.escape(name)}</h3><img alt='{name}' "
                     f"src='data:image/png;base64,{b64}'>")
    page = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><title>"
            f"{html.escape(org)} / {html.escape(ab)}</title><style>{CSS}</style></head><body>"
            + "\n".join(parts) + "</body></html>")
    out.mkdir(parents=True, exist_ok=True)
    (out / "report.html").write_text(page, encoding="utf-8")
    return out / "report.html"


# ---- numbers ---------------------------------------------------------------------------
def numbers(t: Tables, organisms: list[str]) -> pd.DataFrame:
    rows: list[dict] = []

    def put(key, value, table, column, description):
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return
        if isinstance(value, float):
            value = round(value, 4)
        rows.append({"key": key, "value": value, "source_table": table, "column": column,
                     "description": description})

    dec = t.get("panel_decisions")
    if dec is not None:
        put("panel.n_pairs", len(dec), "panel_decisions", "decision", "Organism–antibiotic pairs assessed")
        for d in ("included", "excluded", "blocked"):
            put(f"panel.n_{d}", int((dec["decision"] == d).sum()), "panel_decisions", "decision",
                f"Pairs {d}")
    for org in organisms:
        long = t.get("phenotypes_long", org)
        if long is not None:
            put(f"{org}.n_genomes_phenotyped", int(long["genome_id"].nunique()), "phenotypes_long",
                "genome_id", "Genomes with a phenotype")
            put(f"{org}.n_antibiotics", int(long["antibiotic"].nunique()), "phenotypes_long",
                "antibiotic", "Antibiotics with a phenotype")
        qc = t.get("genome_qc", org)
        if qc is not None:
            put(f"{org}.n_genomes_qc_pass", int(F._bool(qc["pass_overall"]).sum()), "genome_qc",
                "pass_overall", "Genomes that pass QC")
        cl = t.get("lineage_clusters", org)
        if cl is not None:
            put(f"{org}.n_lineages", int(cl["Cluster"].nunique()), "lineage_clusters", "Cluster",
                "PopPUNK lineages")
    included = (dec[(dec["decision"] == panel.INCLUDED) & dec["organism"].isin(organisms)]
                if dec is not None else pd.DataFrame(columns=["organism", "antibiotic"]))
    for org, ab in zip(included["organism"], included["antibiotic"], strict=True):
        mid = f"{org}__{ab}"
        mg = t.get("model_genomes", org, ab)
        if mg is not None:
            put(f"{mid}.n_genomes", len(mg), "model_genomes", "Genome ID", "Genomes of the model")
            put(f"{mid}.n_resistant", int((mg["label"] == 1).sum()), "model_genomes", "label",
                "Resistant genomes")
        rm = t.get("repeat_metrics", org, ab)
        if rm is not None:
            for arm in ("lineage_aware", "lineage_blind"):
                r = rm[rm["arm"] == arm]
                for col in ("roc_auc", "pr_auc", "balanced_accuracy", "mcc"):
                    if col in r and len(r):
                        put(f"{mid}.{col}.{arm}", float(r[col].mean()), "repeat_metrics", col,
                            f"Mean over repeats, {arm.replace('_', '-')} arm")
        grades = t.get("grades_patterns", org, ab)
        if grades is not None:
            put(f"{mid}.n_candidates", len(grades), "grades_patterns", "pattern_id",
                "Candidate patterns")
            for mode in ("allele_aware", "homolog_only"):
                vc = grades[f"{mode}_grade"].value_counts()
                for g in F.GRADES:
                    put(f"{mid}.n_{g}.{mode}", int(vc.get(g, 0)), "grades_patterns",
                        f"{mode}_grade", f"Patterns graded {g}, {mode.replace('_', '-')} mode")
    lp = t.get("label_permutation_models")
    if lp is not None and "q" in lp:
        put("label_permutation.n_significant", int((lp["q"] < float(
            t.config["label_permutation"]["alpha"])).sum()), "label_permutation_models", "q",
            "Models whose label permutation q is below α")
    return pd.DataFrame(rows, columns=["key", "value", "source_table", "column", "description"])


def main():
    ap = argparse.ArgumentParser(description="Figures, model reports and the numbers file.")
    ap.add_argument("--entry", choices=["main", "DOWNLOAD", "CONTEXT", "KB", "CABBAGE"], required=True)
    ap.add_argument("--organisms", required=True)
    args = ap.parse_args()
    config = load_config()
    organisms = args.organisms.split(",")
    t = Tables(config)
    out = resolve_path("reports_dir", config=config)
    n_figs = len(global_figures(t, out))
    for org in organisms:
        n_figs += len(organism_figures(t, org, out / org))
    n_reports = 0
    dec = t.get("panel_decisions")
    if args.entry == "main" and dec is not None:
        inc = dec[(dec["decision"] == panel.INCLUDED) & dec["organism"].isin(organisms)]
        for org, ab in zip(inc["organism"], inc["antibiotic"], strict=True):
            mdir = resolve_path("model_reports_dir", organism=org, antibiotic=ab, config=config)
            figs = model_figures(t, org, ab, mdir)
            n_figs += len(figs)
            model_report(t, org, ab, figs, mdir)
            n_reports += 1
    nums = numbers(t, organisms)
    out.mkdir(parents=True, exist_ok=True)
    nums.to_csv(out / "tez_sayilari.csv", index=False)
    (out / "reports_summary.json").write_text(json.dumps({
        "entry": args.entry, "figures": n_figs, "model_reports": n_reports,
        "numbers": len(nums), "skipped": t.skipped}, indent=2) + "\n")
    print(f"REPORTS — {args.entry}: {n_figs} figures, {n_reports} model reports, "
          f"{len(nums)} numbers")


if __name__ == "__main__":
    main()

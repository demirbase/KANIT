#!/usr/bin/env python3
"""Figures, model reports and the numbers file (plan F4.2), from the steps' CSV tables.

  reports.py --entry main --organisms ecoli,kpneumoniae

At the end of every entry it draws the QC figure of every step whose tables exist
(lib/figures.py; Appendix C), writes for every panel model a self-contained HTML
report (data, cross-validation, label permutation, evidence layers, grades, the best
candidates, apart from them the patterns nominated by association (protocol §14 item 6),
the comparison with genotype-based tools and the model's figures), and
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
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import card_tiers, concordance, contract, panel, registry  # noqa: E402
from lib import figures as F  # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402
from lib.run_metadata import git_commit_hash  # noqa: E402

GRADE_ORDER = {g: i for i, g in enumerate(F.GRADES)}


def _gram_stain(organism_id):
    """The organism's Gram stain from the registry, or None when it has none."""
    try:
        return registry.get_organism(organism_id).get("gram_stain")
    except KeyError:
        return None


def _four_significant(x):
    """Round to four significant digits. A p-value can be far below the fourth
    decimal, where the usual rounding of the numbers file would report it as 0."""
    if x is None or x == 0:
        return x
    return round(x, -int(math.floor(math.log10(abs(x)))) + 3)


class Tables:
    """The contract's tables of one run, read when asked for (None when absent)."""

    def __init__(self, config: dict):
        self.config = config
        self.skipped: list[str] = []

    def get(self, tid: str, organism=None, antibiotic=None, **read_csv) -> pd.DataFrame | None:
        try:
            p = contract.table_path(tid, self.config, organism, antibiotic)
        except KeyError:
            return None
        if not p.exists():
            return None
        return pd.read_csv(p, keep_default_na=True, low_memory=False, **read_csv)

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
    assoc = t.get("association_grades_patterns", org, ab)
    ranks = t.get("association_patterns", org, ab)
    assoc_counts = assoc_best = None
    if assoc is not None and ranks is not None:
        assoc_counts = pd.DataFrame({m: assoc[f"{m}_grade"].value_counts().reindex(
            F.GRADES, fill_value=0) for m in ("allele_aware", "homolog_only")}).reset_index(
            names="grade")
        assoc_best = (assoc.merge(ranks, on="pattern_id", how="left")
                      .assign(_o=lambda d: d["grade"].map(GRADE_ORDER))
                      .sort_values(["_o", "rank"]).head(25)
                      [["pattern_id", "rank", "lrt_pvalue", "n_members", "layers_passed",
                        "allele_aware_card_state", "allele_aware_grade", "homolog_only_grade"]])
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
             "<h2>Patterns nominated by association (secondary analysis)</h2>",
             "<p class='meta'>pyseer-significant patterns that are not candidates (protocol "
             "§14 item 6). They pass the pyseer layer by the way they are chosen, so their "
             "grades are not comparable with those of the candidates.</p>",
             _table(assoc_counts), _table(assoc_best),
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
# Every aggregate the evidence sections of the grading article assert. They were written as
# literals because no key produced them; a number a reader cannot regenerate is not a result.
# `ev` is one row per candidate pattern of every model: the grade tables joined to the
# candidate tables, plus the model id, its organism and its drug class.
def _evidence_numbers(ev: pd.DataFrame, put) -> None:
    layers = ("prevalence", "mda", "cpss", "pyseer")
    for layer in layers:
        ev[layer] = F._bool(ev[layer])
    ev["n_layers_passed"] = ev[list(layers)].sum(axis=1)
    n = len(ev)

    # the layers, one by one and together
    for layer in layers:
        k = int(ev[layer].sum())
        put(f"all.n_pass_{layer}", k, "grades_patterns", layer,
            f"Candidates passing the {layer} layer")
        put(f"all.share_pass_{layer}", k / n, "grades_patterns", layer,
            f"Share (0-1) of the candidates passing the {layer} layer")
    for k in range(len(layers) + 1):
        put(f"all.n_layers_{k}", int((ev["n_layers_passed"] == k).sum()), "grades_patterns",
            "passes", f"Candidates passing exactly {k} statistical layer(s)")
    combos = ev.groupby([ev[l] for l in layers]).size().sort_values(ascending=False)
    for rank, (flags, size) in enumerate(combos.items(), start=1):
        passed = [l for l, f in zip(layers, flags, strict=True) if f]
        if rank > 3 or not passed:
            continue
        put(f"all.n_layer_combo_{'_and_'.join(passed)}", int(size), "grades_patterns", "passes",
            f"Candidates whose passing layers are exactly {' + '.join(passed)}")

    # the two grading modes: which patterns move, and where to
    moved = ev[ev["allele_aware_grade"] != ev["homolog_only_grade"]]
    put("all.n_grade_differs_between_modes", len(moved), "grades_patterns", "homolog_only_grade",
        "Candidates graded differently in the allele-aware and homolog-only modes")
    put("all.share_grade_differs_between_modes", len(moved) / n, "grades_patterns",
        "homolog_only_grade", "Share (0-1) of the candidates whose grade differs between modes")
    for (a, b), size in moved.groupby(["allele_aware_grade", "homolog_only_grade"]).size().items():
        put(f"all.n_mode_change_{a}_to_{b}", int(size), "grades_patterns", "homolog_only_grade",
            f"Candidates graded {a} allele-aware and {b} homolog-only")

    # what the gain ranking alone reports (the top-gain patterns of every model)
    top = ev[ev["gain_rank"].notna()]
    if len(top):
        nt = len(top)
        put("all.n_top_gain_patterns", nt, "candidates", "gain_rank",
            "Patterns in the top total gain of a final model, summed over models")
        states = {"b": "b", "card_hit_without_b": "card_hit_without_b", "no_card_hit": "no_card_hit"}
        for name, state in states.items():
            k = int((top["allele_aware_card_state"] == state).sum())
            put(f"all.top_gain_n_{name}", k, "grades_patterns", "allele_aware_card_state",
                f"Top-gain patterns whose CARD state is {state}")
            put(f"all.top_gain_share_{name}", k / nt, "grades_patterns",
                "allele_aware_card_state", f"Share (0-1) of the top-gain patterns with {state}")
        first = top.loc[top.groupby("model")["gain_rank"].idxmin()]
        put("all.n_models_top_feature_b",
            int((first["allele_aware_card_state"] == "b").sum()), "grades_patterns",
            "allele_aware_card_state", "Models whose single highest-gain pattern has b = 1")
        hit = top[top["allele_aware_card_state"] == "card_hit_without_b"]
        for reason, size in hit["allele_aware_card_reasons"].value_counts().items():
            put(f"all.top_gain_n_reason_{reason}", int(size), "grades_patterns",
                "allele_aware_card_reasons",
                f"Top-gain patterns denied b = 1 because of {reason}")
            put(f"all.top_gain_share_reason_{reason}", int(size) / len(hit), "grades_patterns",
                "allele_aware_card_reasons",
                f"Share (0-1) of those patterns denied b = 1 because of {reason}")
        blank = top[top["allele_aware_card_state"] == "no_card_hit"]
        failed = blank[~blank["pyseer"]]
        put("all.top_gain_no_hit_n_pyseer_fail", len(failed), "grades_patterns", "pyseer",
            "Top-gain patterns without a CARD hit that do not pass the mixed model")
        if len(blank):
            put("all.top_gain_no_hit_share_pyseer_fail", len(failed) / len(blank),
                "grades_patterns", "pyseer",
                "Share (0-1) of the top-gain patterns without a CARD hit failing the mixed model")
        put("all.top_gain_no_hit_pyseer_fail_n_prevalence_pass",
            int(failed["prevalence"].sum()), "grades_patterns", "prevalence",
            "Of those, the ones that still pass the prevalence difference")
        for g in F.GRADES:
            k = int((top["grade"] == g).sum())
            put(f"all.top_gain_share_{g}", k / nt, "grades_patterns", "grade",
                f"Share (0-1) of the top-gain patterns graded {g}")
            put(f"all.n_{g}_in_top_gain", k, "grades_patterns", "grade",
                f"Patterns graded {g} that a model ranks in its top total gain")
            put(f"all.n_{g}_outside_top_gain", int((ev["grade"] == g).sum()) - k,
                "grades_patterns", "grade",
                f"Patterns graded {g} reaching the candidate set only by stability selection")

    # the grade a candidate reaches without a CARD hit, with one, or with b = 1
    cand = ev[ev["grade"] == "candidate"]
    if len(cand):
        for name, state in (("no_card_hit", "no_card_hit"),
                            ("card_hit_without_b", "card_hit_without_b"), ("b", "b")):
            put(f"all.n_candidate_grade_{name}",
                int((cand["allele_aware_card_state"] == state).sum()), "grades_patterns",
                "allele_aware_card_state", f"Candidate-graded patterns whose CARD state is {state}")
        hit = cand[cand["allele_aware_card_state"] == "card_hit_without_b"]
        for reason, size in hit["allele_aware_card_reasons"].value_counts().items():
            put(f"all.n_candidate_grade_reason_{reason}", int(size), "grades_patterns",
                "allele_aware_card_reasons", f"Candidate-graded patterns denied b = 1 by {reason}")

    novel = ev[ev["grade"] == "strong_novel"]
    if len(novel):
        put("all.n_models_with_strong_novel", int(novel["model"].nunique()), "grades_patterns",
            "grade", "Models contributing at least one strong_novel pattern")
        for k in (3, 4):
            put(f"all.n_strong_novel_layers_{k}", int((novel["n_layers_passed"] == k).sum()),
                "grades_patterns", "passes", f"strong_novel patterns passing {k} layers")

    # H2: the share of every model's stable patterns that reach b = 1
    stable = ev[F._bool(ev["cpss_stable"])] if "cpss_stable" in ev else ev.iloc[:0]
    if len(stable):
        put("all.n_stable_patterns", len(stable), "candidates", "cpss_stable",
            "Patterns the stability selection calls stable, summed over models")
        k = int((stable["allele_aware_card_state"] == "b").sum())
        put("all.n_stable_b", k, "grades_patterns", "allele_aware_card_state",
            "Of those, the ones matching a curated determinant of the modelled drug class")
        put("all.stable_b_share", k / len(stable), "grades_patterns", "allele_aware_card_state",
            "Share (0-1) of the stable patterns with b = 1, pooled over models")
        for mode in ("allele_aware", "homolog_only"):
            share = (stable.assign(hit=stable[f"{mode}_card_state"] == "b")
                     .groupby("model")["hit"].mean())
            put(f"all.stable_b_median_share.{mode}", float(share.median()), "grades_patterns",
                f"{mode}_card_state",
                f"Median over models of the share of stable patterns with b = 1, {mode} mode")
            if mode == "allele_aware":
                put("all.stable_b_share_q1", float(share.quantile(0.25)), "grades_patterns",
                    "allele_aware_card_state", "First quartile of that share over models")
                put("all.stable_b_share_q3", float(share.quantile(0.75)), "grades_patterns",
                    "allele_aware_card_state", "Third quartile of that share over models")
                put("all.n_models_with_stable_patterns", int(share.size), "candidates",
                    "cpss_stable", "Models with at least one stable pattern")
                put("all.n_models_no_stable_b", int((share == 0).sum()), "grades_patterns",
                    "allele_aware_card_state",
                    "Models where no stable pattern matches a curated determinant")

    # grade shares per organism and per drug class, which 3.3 reads off the figure
    for col, label in (("organism", "organism"), ("drug_class", "drug class")):
        if col not in ev:
            continue
        for key, part in ev.groupby(col):
            slug = str(key).replace(" ", "_").replace("/", "_")
            put(f"{slug}.n_candidates_in_group", len(part), "grades_patterns", "pattern_id",
                f"Candidate patterns of this {label}")
            for g in F.GRADES:
                put(f"{slug}.n_{g}", int((part["grade"] == g).sum()), "grades_patterns", "grade",
                    f"Candidate patterns of this {label} graded {g}")
                put(f"{slug}.share_{g}", float((part["grade"] == g).mean()), "grades_patterns",
                    "grade", f"Share (0-1) of this {label}'s candidates graded {g}")


# hypotheses.json (protocol §10) -> keys. A key names what its test measures, not the
# test's number: the protocol's H3 and H7 are the article's H4 and H3, so a number in the
# key would be right for one text and wrong for the other.
# (summary, field, key, table, column, description)
HYPOTHESIS_KEYS = [
    ("H1", "n_models", "stable.n_models", "h1_models", "model_id", "H1: models tested"),
    ("H1", "n_with_stable", "stable.n_with_stable", "h1_models", "has_stable",
     "H1: models with at least one stable pattern"),
    ("H1", "supported", "stable.supported", "h1_models", "has_stable",
     "H1: 1 if every model has a stable pattern, else 0"),
    ("H2", "n_defined", "curated.n_models", "h2_models", "share_b_allele_aware",
     "H2: models with stable patterns, the ones the median is taken over"),
    ("H2", "median", "curated.median_share", "h2_models", "share_b_allele_aware",
     "H2: median over models of the share of stable patterns with b = 1"),
    ("H2", "median_homolog_only", "curated.median_share_homolog_only", "h2_models",
     "share_b_homolog_only", "H2: that median in the homolog-only mode"),
    ("H2", "supported", "curated.supported", "h2_models", "share_b_allele_aware",
     "H2: 1 if the median is at least 0.40, else 0"),
    ("H3", "n_within", "class_overlap.n_within", "h3_pairs", "same_class",
     "Protocol H3: within-class model pairs with a defined overlap"),
    ("H3", "n_cross", "class_overlap.n_cross", "h3_pairs", "same_class",
     "Protocol H3: cross-class model pairs with a defined overlap"),
    ("H3", "median_within", "class_overlap.median_within", "h3_pairs", "overlap",
     "Protocol H3: median overlap coefficient of the within-class pairs"),
    ("H3", "median_cross", "class_overlap.median_cross", "h3_pairs", "overlap",
     "Protocol H3: median overlap coefficient of the cross-class pairs"),
    ("H3", "p", "class_overlap.p", "h3_pairs", "overlap",
     "Protocol H3: one-sided Mann-Whitney U p, within-class above cross-class"),
    ("H3", "supported", "class_overlap.supported", "h3_pairs", "overlap",
     "Protocol H3: 1 if p < 0.05, else 0"),
    ("H6", "n_acquired", "acquired_vs_mutational.n_acquired", "h6_models", "class",
     "H6 (explanatory): acquired phenotypes with an H2 share"),
    ("H6", "n_mutational", "acquired_vs_mutational.n_mutational", "h6_models", "class",
     "H6 (explanatory): mutational phenotypes with an H2 share"),
    ("H6", "n_mixed", "acquired_vs_mutational.n_mixed", "h6_models", "class",
     "H6 (explanatory): mixed phenotypes"),
    ("H6", "p_acquired_gt_mutational", "acquired_vs_mutational.p", "h6_models", "h2_share",
     "H6 (explanatory): one-sided Mann-Whitney U p, acquired above mutational"),
]
# protocol H7 at each threshold of expected prevalence among the resistant genomes
EXPECTED_THRESHOLDS = {"0.1": "thr10", "0.05": "thr5", "0.2": "thr20"}
EXPECTED_KEYS = [
    ("n_expected", "n_expected", "expected determinants"),
    ("n_recovered", "n_recovered", "expected determinants a candidate maps to"),
    ("share_confirmed", "share_confirmed", "share (0-1) of those whose best grade is confirmed"),
    ("supported_primary", "supported_primary", "1 if that share is at least 0.5, else 0"),
    ("secondary_p", "secondary_p",
     "one-sided Fisher p, confirmed more often among candidates of an expected determinant"),
    ("supported_secondary", "supported_secondary", "1 if that p is below 0.05, else 0"),
]
SECONDARY_CELLS = [((0, 0), "expected_confirmed", "candidates of an expected determinant, confirmed"),
                   ((0, 1), "expected_not_confirmed",
                    "candidates of an expected determinant, not confirmed"),
                   ((1, 0), "other_confirmed",
                    "candidates with a CARD hit to another family, confirmed"),
                   ((1, 1), "other_not_confirmed",
                    "candidates with a CARD hit to another family, not confirmed")]


def _hypothesis_numbers(summary: dict, put, pre: str, label: str = "") -> None:
    """The verdicts and statistics of one set of hypothesis summaries under pre."""
    def val(x):
        return int(x) if isinstance(x, bool) else x

    for test, field, key, table, column, text in HYPOTHESIS_KEYS:
        x = val((summary.get(test) or {}).get(field))
        if field.startswith("p"):
            put(f"{pre}.{key}", _four_significant(x), table, column, label + text, digits=None)
        else:
            put(f"{pre}.{key}", x, table, column, label + text)
    for thr, tag in EXPECTED_THRESHOLDS.items():
        s = (summary.get("H7") or {}).get(thr) or {}
        for field, key, text in EXPECTED_KEYS:
            x = val(s.get(field))
            p_value = field == "secondary_p"
            put(f"{pre}.expected.{key}.{tag}", _four_significant(x) if p_value else x,
                "h7_determinants", "best_grade", f"{label}Protocol H7, threshold {thr}: {text}",
                digits=None if p_value else 6 if field == "share_confirmed" else "auto")
        table = s.get("secondary_table")
        if table:
            for (i, j), key, text in SECONDARY_CELLS:
                put(f"{pre}.expected.secondary_{key}.{tag}", int(table[i][j]), "h7_determinants",
                    "best_grade", f"{label}Protocol H7, threshold {thr}: {text}")


def _flagged_models(lp: pd.DataFrame, models: list[str], alpha: float) -> list[str] | None:
    """Models whose label permutation q is at or above α (§8.6); None unless every model
    has its row, since the q values are corrected across all of them."""
    if "q" not in lp:
        return None
    q = dict(zip(lp["organism"] + "__" + lp["antibiotic"], lp["q"], strict=True))
    if any(m not in q for m in models):
        return None
    return sorted(m for m in models if pd.notna(q[m]) and q[m] >= alpha)


def _permutation_sensitivity(ev: pd.DataFrame, flagged: list[str], put) -> None:
    """Protocol §14 item 9: the grade distribution without the flagged models, and the
    layers each flagged model's candidates pass (an illustration, not a test). The grades
    of a flagged model are its {model}.n_{grade} keys. The hypotheses are recomputed in
    hypotheses.json."""
    layers = ("prevalence", "mda", "cpss", "pyseer")
    rest = ev[~ev["model"].isin(flagged)]
    if flagged and len(rest):
        pre = "lp_sensitivity"
        put(f"{pre}.n_models", int(rest["model"].nunique()), "grades_patterns", "pattern_id",
            "Models left once the models flagged by label permutation are excluded")
        put(f"{pre}.n_candidate_patterns", len(rest), "grades_patterns", "pattern_id",
            "Candidate patterns of those models")
        for g in F.GRADES:
            k = int((rest["grade"] == g).sum())
            put(f"{pre}.n_{g}", k, "grades_patterns", "grade",
                f"Candidate patterns graded {g} without the flagged models, reported rule")
            put(f"{pre}.share_{g}", k / len(rest), "grades_patterns", "grade",
                f"Share (0-1) of those candidate patterns graded {g}, reported rule")
        put(f"{pre}.weak_or_none_share", float(rest["grade"].isin(("weak", "none")).mean()),
            "grades_patterns", "grade",
            "Share (0-1) of those candidate patterns graded weak or none, reported rule")
        for layer in layers:
            put(f"{pre}.share_pass_{layer}", float(F._bool(rest[layer]).mean()),
                "grades_patterns", layer,
                f"Share (0-1) of those candidate patterns passing the {layer} layer")
    for mid in flagged:
        part = ev[ev["model"] == mid]
        for layer in layers:
            put(f"{mid}.n_pass_{layer}", int(F._bool(part[layer]).sum()), "grades_patterns",
                layer, f"Candidates of this flagged model passing the {layer} layer")


def card_tier_table(t: Tables, models: list[tuple[str, str]]) -> pd.DataFrame:
    """The candidates with a CARD hit of every model, with the RGI tier of their hits
    (protocol §14 item 10). Empty when no organism has its hit table."""
    parts = []
    for org in dict.fromkeys(o for o, _ in models):
        hits = t.get("rgi_hits", org, usecols=["genome_id", "aro", "cut_off"],
                     dtype={"genome_id": str, "aro": str})
        if hits is None:
            continue
        tiers = card_tiers.pair_tiers(hits)
        for o, ab in models:
            if o != org:
                continue
            pats = t.get("grades_patterns", org, ab)
            units = t.get("grades_unitigs", org, ab, dtype={"located_genomes": str,
                          f"{card_tiers.RULE}_card_aros": str})
            if pats is None or units is None:
                continue
            parts.append(card_tiers.model_table(pats, units, tiers)
                         .assign(model=f"{org}__{ab}", organism=org))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _card_tier_numbers(table: pd.DataFrame, put) -> None:
    src = "grades_unitigs"
    tiered = table[table["tier"].notna() & table["group"].notna()]
    put("card_tier.n_patterns", len(tiered), src, "located_genomes",
        "Candidates with a CARD hit and an RGI tier (§14 item 10)")
    put("card_tier.n_patterns_without_tier", int(table["tier"].isna().sum()), src,
        "located_genomes", "Candidates with a CARD hit whose every (genome, ARO) pair is "
        "ambiguous or has no hit")
    put("card_tier.n_models", int(tiered["model"].nunique()), src, "pattern_id",
        "Models with at least one such candidate")
    aros = {a for x in tiered["aros"].dropna() for a in str(x).split(";")}
    put("card_tier.n_aros", len(aros), "rgi_hits", "aro", "AROs behind those candidates")
    for c, what in (("n_perfect", "Perfect"), ("n_strict", "Strict"),
                    ("n_ambiguous", "ambiguous (hits of the ARO in both tiers)"),
                    ("n_no_hit", "without a Perfect or Strict hit of the ARO")):
        put(f"card_tier.n_pairs_{c[2:]}", int(table[c].fillna(0).sum()), "rgi_hits", "cut_off",
            f"(genome, ARO) pairs of the candidates, {what}")
    for g in card_tiers.GROUPS:
        s = tiered[tiered["group"] == g]
        for tier in card_tiers.TIER_ORDER:
            put(f"card_tier.n_{g}_{tier.lower()}", int((s["tier"] == tier).sum()), "rgi_hits",
                "cut_off", f"Candidates of group {g} whose tier is {tier}")
        if len(s):
            put(f"card_tier.share_perfect_{g}", float((s["tier"] == "Perfect").mean()),
                "rgi_hits", "cut_off", f"Share (0-1) of the group {g} whose tier is Perfect")
    for org, s in tiered.groupby("organism"):
        for g in card_tiers.COMPARED:
            sg = s[s["group"] == g]
            put(f"card_tier.{org}.n_{g}", len(sg), "rgi_hits", "cut_off",
                f"Candidates of group {g} with a tier")
            if len(sg):
                put(f"card_tier.{org}.share_perfect_{g}", float((sg["tier"] == "Perfect").mean()),
                    "rgi_hits", "cut_off", f"Share (0-1) of the group {g} whose tier is Perfect")
    c = card_tiers.compare(tiered)
    put("card_tier.compare.n_confirmed", c["n_confirmed"], src, "allele_aware_grade",
        "confirmed candidates with b = 1 and a tier, in the comparison")
    put("card_tier.compare.n_candidate", c["n_candidate_b"], src, "allele_aware_grade",
        "candidate candidates with b = 1 and a tier, in the comparison")
    put("card_tier.compare.n_organisms", c["n_strata"], "rgi_hits", "cut_off",
        "Organisms with both grades and both outcomes, which the test is computed over")
    put("card_tier.compare.odds_ratio", c["odds_ratio"], "rgi_hits", "cut_off",
        "Mantel–Haenszel odds ratio of a Perfect tier, confirmed vs candidate, by organism")
    put("card_tier.compare.ci_low", c["ci_low"], "rgi_hits", "cut_off",
        "Lower 95% bound of that odds ratio (Robins–Breslow–Greenland)")
    put("card_tier.compare.ci_high", c["ci_high"], "rgi_hits", "cut_off",
        "Upper 95% bound of that odds ratio (Robins–Breslow–Greenland)")
    if c["p"] is not None:
        put("card_tier.compare.p", _four_significant(c["p"]), "rgi_hits", "cut_off",
            "Cochran–Mantel–Haenszel p, two-sided, no continuity correction", digits=None)


def numbers(t: Tables, organisms: list[str]) -> pd.DataFrame:
    rows: list[dict] = []

    def put(key, value, table, column, description, digits="auto"):
        # A share rounded to four decimals cannot be rendered as a percentage with one
        # decimal: 0.074471 becomes 0.0745, which reads as 7.5 % instead of 7.4 %.
        # digits=None keeps the value as it is: a p-value, already rounded to four
        # significant digits, would read 0 at four decimals.
        if digits == "auto":
            digits = 6 if description.startswith("Share (0-1)") else 4
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return
        if isinstance(value, float) and digits is not None:
            value = round(value, digits)
        rows.append({"key": key, "value": value, "source_table": table, "column": column,
                     "description": description})

    min_support = int((t.config.get("unitig") or {}).get("min_support") or 10)
    dec = t.get("panel_decisions")
    if dec is not None:
        put("panel.n_pairs", len(dec), "panel_decisions", "decision", "Organism–antibiotic pairs assessed")
        for d in ("included", "excluded", "blocked"):
            put(f"panel.n_{d}", int((dec["decision"] == d).sum()), "panel_decisions", "decision",
                f"Pairs {d}")
        # 3.1 splits the excluded pairs by why they fell out, and quotes the minority class of
        # the pairs that stayed; both were literals because no key carried them.
        exc = dec[dec["decision"] == "excluded"]
        if len(exc) and "reason" in exc:
            short = exc["reason"].str.startswith("minority")
            put("panel.n_excluded_small_minority", int(short.sum()), "panel_decisions", "reason",
                "Pairs excluded because the minority class is below the threshold")
            put("panel.n_excluded_other_reason", int((~short).sum()), "panel_decisions", "reason",
                "Pairs excluded for any other reason, e.g. a label naming no single drug")
        # The unitig/k-mer article weighs a relative frequency floor against min_support. A
        # floor of q is stricter exactly where q * n_eligible > min_support, so the number of
        # pairs it would bind is what that argument turns on.
        inc_n = pd.to_numeric(
            dec.loc[dec["decision"] == panel.INCLUDED, "n_eligible"], errors="coerce").dropna() \
            if "n_eligible" in dec else pd.Series(dtype=float)
        if len(inc_n):
            for name, q in (("1pct", 0.01), ("5pct", 0.05)):
                put(f"panel.n_pairs_relative_floor_stricter_{name}",
                    int((inc_n * q > min_support).sum()), "panel_decisions", "n_eligible",
                    f"Included pairs where a {q:.0%} relative floor binds tighter than "
                    f"min_support = {min_support}")
            put("panel.n_eligible_median", float(inc_n.median()), "panel_decisions", "n_eligible",
                "Median genomes of the included pairs after QC and lineage assignment")
            put("panel.n_eligible_min", int(inc_n.min()), "panel_decisions", "n_eligible",
                "Fewest genomes of any included pair")
            put("panel.n_eligible_max", int(inc_n.max()), "panel_decisions", "n_eligible",
                "Most genomes of any included pair")
        inc_min = pd.to_numeric(
            dec.loc[dec["decision"] == panel.INCLUDED, "minority"], errors="coerce").dropna()
        if len(inc_min):
            put("panel.median_minority_included", float(inc_min.median()), "panel_decisions",
                "minority", "Median minority class of the included pairs")
            put("panel.min_minority_included", int(inc_min.min()), "panel_decisions", "minority",
                "Smallest minority class of the included pairs")
            put("panel.max_minority_included", int(inc_min.max()), "panel_decisions", "minority",
                "Largest minority class of the included pairs")
    comparisons: list[pd.DataFrame] = []
    genomes_qc_pass = 0
    for org in organisms:
        long = t.get("phenotypes_long", org)
        if long is not None:
            put(f"{org}.n_genomes_phenotyped", int(long["genome_id"].nunique()), "phenotypes_long",
                "genome_id", "Genomes with a phenotype")
            put(f"{org}.n_antibiotics", int(long["antibiotic"].nunique()), "phenotypes_long",
                "antibiotic", "Antibiotics with a phenotype")
        qc = t.get("genome_qc", org)
        if qc is not None:
            n_pass = int(F._bool(qc["pass_overall"]).sum())
            put(f"{org}.n_genomes_qc_pass", n_pass, "genome_qc",
                "pass_overall", "Genomes that pass QC")
            genomes_qc_pass += n_pass
        if dec is not None:
            mine = dec[dec["organism"] == org]
            if len(mine):
                put(f"{org}.n_pairs_assessed", len(mine), "panel_decisions", "decision",
                    "Pairs of this organism the panel rule assessed")
                put(f"{org}.n_pairs_included", int((mine["decision"] == panel.INCLUDED).sum()),
                    "panel_decisions", "decision", "Pairs of this organism that entered the panel")
                if "minority" in mine:
                    put(f"{org}.max_minority", int(mine["minority"].max()), "panel_decisions",
                        "minority", "Largest minority class among this organism's pairs")
                for r in mine.itertuples(index=False):
                    for col, what in (("n_eligible", "Genomes after QC and lineage assignment"),
                                      ("n_tested", "Genomes with a phenotype for this pair")):
                        v = pd.to_numeric(getattr(r, col, None), errors="coerce")
                        if pd.notna(v):
                            put(f"{r.organism}__{r.antibiotic}.{col}", int(v), "panel_decisions",
                                col, what)
        cl = t.get("lineage_clusters", org)
        if cl is not None:
            put(f"{org}.n_lineages", int(cl["Cluster"].nunique()), "lineage_clusters", "Cluster",
                "PopPUNK lineages")
        comp = t.get("external_comparison", org)
        if comp is not None:
            comparisons.append(comp)
    included = (dec[(dec["decision"] == panel.INCLUDED) & dec["organism"].isin(organisms)]
                if dec is not None else pd.DataFrame(columns=["organism", "antibiotic"]))
    if dec is not None:
        put("panel.n_organisms_assessed", int(dec["organism"].nunique()),
            "panel_decisions", "organism", "Organisms the panel rule assessed")
        put("panel.n_organisms_with_models", int(included["organism"].nunique()),
            "panel_decisions", "organism", "Organisms with at least one included pair")
        put("panel.n_antibiotics_included", int(included["antibiotic"].nunique()),
            "panel_decisions", "antibiotic", "Antibiotics of the included pairs")
        if "drug_class" in included:
            put("panel.n_drug_classes_included", int(included["drug_class"].nunique()),
                "panel_decisions", "drug_class", "Drug classes of the included pairs")
    genomes_in_models: set[str] = set()
    evidence: list[pd.DataFrame] = []
    drug_class = (dict(zip(included["organism"] + "__" + included["antibiotic"],
                           included["drug_class"], strict=True))
                  if "drug_class" in included else {})
    n_candidate_patterns = 0
    n_graded_primary = 0
    n_weak_or_none = 0
    grade_counts: dict[str, int] = {g: 0 for g in F.GRADES}
    for org, ab in zip(included["organism"], included["antibiotic"], strict=True):
        mid = f"{org}__{ab}"
        mg = t.get("model_genomes", org, ab)
        if mg is not None:
            put(f"{mid}.n_genomes", len(mg), "model_genomes", "Genome ID", "Genomes of the model")
            genomes_in_models.update(mg["Genome ID"].astype(str))
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
            n_candidate_patterns += len(grades)
            if "grade" in grades:
                n_graded_primary += len(grades)
                n_weak_or_none += int(grades["grade"].isin(("weak", "none")).sum())
                vc_primary = grades["grade"].value_counts()
                for g in F.GRADES:
                    grade_counts[g] += int(vc_primary.get(g, 0))
            for mode in ("allele_aware", "homolog_only"):
                vc = grades[f"{mode}_grade"].value_counts()
                for g in F.GRADES:
                    put(f"{mid}.n_{g}.{mode}", int(vc.get(g, 0)), "grades_patterns",
                        f"{mode}_grade", f"Patterns graded {g}, {mode.replace('_', '-')} mode")
        cpss = t.get("cpss_patterns", org, ab)
        if cpss is not None and "pi" in cpss and len(cpss):
            pi = pd.to_numeric(cpss["pi"], errors="coerce")
            put(f"{mid}.max_cpss_pi", float(pi.max()), "cpss_patterns", "pi",
                "Highest selection frequency the stability selection reached")
            put(f"{mid}.n_patterns_cpss_pi_positive", int((pi > 0).sum()), "cpss_patterns", "pi",
                "Patterns the stability selection picked at least once")
        cands = t.get("candidates", org, ab)
        if grades is not None and cands is not None and "grade" in grades:
            evidence.append(grades.merge(cands, on="pattern_id", how="left")
                            .assign(model=mid, organism=org,
                                    drug_class=drug_class.get(mid, "unmapped")))
        assoc = t.get("association_grades_patterns", org, ab)
        if assoc is not None:
            put(f"{mid}.association.n_patterns", len(assoc), "association_grades_patterns",
                "pattern_id", "Patterns nominated by association (secondary analysis)")
            vc = assoc["grade"].value_counts()
            for g in F.GRADES:
                put(f"{mid}.association.n_{g}", int(vc.get(g, 0)), "association_grades_patterns",
                    "grade", f"Patterns nominated by association graded {g}, primary mode")
    if evidence:
        ev = pd.concat(evidence, ignore_index=True)
        _evidence_numbers(ev, put)
        lp = t.get("label_permutation_models")
        flagged = (_flagged_models(lp, sorted(ev["model"].unique()),
                                   float(t.config["label_permutation"]["alpha"]))
                   if lp is not None else None)
        if flagged is not None:
            _permutation_sensitivity(ev, flagged, put)
        tier_table = card_tier_table(t, sorted({(m.split("__")[0], m.split("__", 1)[1])
                                                for m in ev["model"]}))
        if len(tier_table):
            _card_tier_numbers(tier_table, put)
    # 3.1 quotes the spread of model performance and of the gap between the two arms. Both are
    # medians over models, so they belong here rather than in any one model's keys.
    auc: dict[str, dict[str, float]] = {"lineage_aware": {}, "lineage_blind": {}}
    for row in rows:
        for arm in auc:
            if row["key"].endswith(f".roc_auc.{arm}"):
                auc[arm][row["key"].rsplit(".", 2)[0]] = row["value"]
    if auc["lineage_aware"]:
        a_vals = pd.Series(auc["lineage_aware"])
        put("models.roc_auc_median.lineage_aware", float(a_vals.median()), "repeat_metrics",
            "roc_auc", "Median over models of the lineage-aware ROC-AUC")
        put("models.roc_auc_min.lineage_aware", float(a_vals.min()), "repeat_metrics", "roc_auc",
            "Lowest lineage-aware ROC-AUC of any model")
        put("models.roc_auc_max.lineage_aware", float(a_vals.max()), "repeat_metrics", "roc_auc",
            "Highest lineage-aware ROC-AUC of any model")
        shared = sorted(set(auc["lineage_aware"]) & set(auc["lineage_blind"]))
        if shared:
            gap = pd.Series({m: auc["lineage_blind"][m] - auc["lineage_aware"][m]
                             for m in shared})
            put("models.arm_gap_median", float(gap.median()), "repeat_metrics", "roc_auc",
                "Median over models of the random-split minus lineage-aware ROC-AUC")
            put("models.arm_gap_max", float(gap.max()), "repeat_metrics", "roc_auc",
                "Largest such gap")
            put("models.arm_gap_min", float(gap.min()), "repeat_metrics", "roc_auc",
                "Smallest such gap; above zero means no model is an exception")
            put("models.n_arm_gap_positive", int((gap > 0).sum()), "repeat_metrics", "roc_auc",
                "Models where the random split gives the higher ROC-AUC")
            by_organism = {}
            for org in organisms:
                mine = gap[[m for m in gap.index if m.startswith(f"{org}__")]]
                if len(mine):
                    by_organism[org] = float(mine.median())
                    put(f"{org}.arm_gap_median", float(mine.median()), "repeat_metrics",
                        "roc_auc", "Median arm gap over this organism's models")
                    put(f"{org}.arm_gap_max", float(mine.max()), "repeat_metrics", "roc_auc",
                        "Largest arm gap among this organism's models")
                own = pd.Series({m: v for m, v in auc["lineage_aware"].items()
                                 if m.startswith(f"{org}__")})
                if len(own):
                    put(f"{org}.roc_auc_median.lineage_aware", float(own.median()),
                        "repeat_metrics", "roc_auc",
                        "Median lineage-aware ROC-AUC of this organism's models")
            if by_organism:
                put("models.arm_gap_median_by_organism_min", min(by_organism.values()),
                    "repeat_metrics", "roc_auc",
                    "Smallest of the per-organism median arm gaps")
                put("models.arm_gap_median_by_organism_max", max(by_organism.values()),
                    "repeat_metrics", "roc_auc", "Largest of the per-organism median arm gaps")
    if genomes_in_models:
        put("all.n_genomes_in_models", len(genomes_in_models), "model_genomes", "Genome ID",
            "Distinct genomes entering at least one model")
    if n_candidate_patterns:
        put("all.n_candidate_patterns", n_candidate_patterns, "grades_patterns", "pattern_id",
            "Candidate patterns, summed over models (pattern numbers are model-specific)")
    if n_graded_primary:
        put("all.weak_or_none_share", n_weak_or_none / n_graded_primary, "grades_patterns",
            "grade", "Share (0-1) of the candidate patterns graded weak or none, reported rule")
        for g in F.GRADES:
            put(f"all.n_{g}", grade_counts[g], "grades_patterns", "grade",
                f"Candidate patterns graded {g}, summed over models, reported rule")
            put(f"all.share_{g}", grade_counts[g] / n_graded_primary, "grades_patterns", "grade",
                f"Share (0-1) of the candidate patterns graded {g}, reported rule")
    if genomes_qc_pass:
        put("all.n_genomes_qc_pass", genomes_qc_pass, "genome_qc", "pass_overall",
            "Genomes that pass QC, summed over organisms")
    if comparisons:
        comp_all = pd.concat(comparisons, ignore_index=True)
        ba = comp_all.pivot_table(index="model_id", columns="tool",
                                  values="balanced_accuracy", aggfunc="first")
        ok = comp_all.pivot_table(index="model_id", columns="tool",
                                  values="assessable", aggfunc="first")
        if "model" in ba and "model" in ok:
            for tool in sorted(c for c in ba.columns if c != "model"):
                both = (ok["model"] == 1) & (ok[tool] == 1)
                d = (ba.loc[both, "model"] - ba.loc[both, tool]).dropna()
                s = concordance.sign_test(d)
                pre = f"compare.{tool}"
                put(f"{pre}_n_assessable", int(both.sum()), "external_comparison", "assessable",
                    f"Pairs both the model and {tool} can assess")
                put(f"{pre}_n_model_better", s["n_positive"], "external_comparison",
                    "balanced_accuracy", f"Pairs where the model beats {tool}")
                put(f"{pre}_n_tool_better", s["n_negative"], "external_comparison",
                    "balanced_accuracy", f"Pairs where {tool} beats the model")
                put(f"{pre}_n_tied", s["n_tied"], "external_comparison", "balanced_accuracy",
                    f"Pairs where the model and {tool} score the same (dropped by the sign test)")
                put(f"{pre}_p", _four_significant(s["p_value"]), "external_comparison",
                    "balanced_accuracy",
                    f"Two-sided sign test of balanced accuracy, model vs {tool}, ties dropped",
                    digits=None)
                # the same test with the ties counted for one side: 3.5 reports how far the
                # result turns on a convention, so both alternatives need a key of their own
                put(f"{pre}_p_tie_to_model", _four_significant(s["p_tie_to_ours"]),
                    "external_comparison", "balanced_accuracy",
                    "That sign test with a tie counted for the model, not dropped", digits=None)
                put(f"{pre}_p_tie_to_tool", _four_significant(s["p_tie_to_theirs"]),
                    "external_comparison", "balanced_accuracy",
                    f"That sign test with a tie counted for {tool}, not dropped", digits=None)
                put(f"{pre}_median_ba_model", float(ba.loc[both, "model"].median()),
                    "external_comparison", "balanced_accuracy",
                    f"Median balanced accuracy of the model, pairs shared with {tool}")
                put(f"{pre}_median_ba_tool", float(ba.loc[both, tool].median()),
                    "external_comparison", "balanced_accuracy",
                    f"Median balanced accuracy of {tool} on those pairs")
                pd_ = concordance.paired_difference(
                    d, [i.split("__")[0] for i in d.index],
                    n_boot=int(t.config["cv"]["n_bootstrap"]))
                put(f"{pre}_median_diff", pd_["median"], "external_comparison",
                    "balanced_accuracy",
                    f"Median paired difference in balanced accuracy, model minus {tool}")
                put(f"{pre}_diff_ci_low", pd_["ci_low"], "external_comparison",
                    "balanced_accuracy",
                    "Lower 95% bound of that difference, organism-clustered bootstrap (§6.3)")
                put(f"{pre}_diff_ci_high", pd_["ci_high"], "external_comparison",
                    "balanced_accuracy",
                    "Upper 95% bound of that difference, organism-clustered bootstrap (§6.3)")
                put(f"{pre}_n_organisms", pd_["n_clusters"], "external_comparison", "model_id",
                    f"Organisms contributing a pair shared with {tool}")
                put(f"{pre}_p_signed_rank",
                    _four_significant(concordance.signed_rank_p(d)), "external_comparison",
                    "balanced_accuracy",
                    f"Signed-rank p, model vs {tool}; sensitivity only, pairs not independent",
                    digits=None)
                rows_t = comp_all[(comp_all["tool"] == tool) & (comp_all["assessable"] == 1)]
                if {"sensitivity", "specificity"} <= set(rows_t.columns):
                    put(f"{pre}_n_calls_all_resistant",
                        int(((rows_t["sensitivity"] == 1) & (rows_t["specificity"] == 0)).sum()),
                        "external_comparison", "sensitivity",
                        f"Models where {tool} calls every genome resistant")
                gram = rows_t["model_id"].str.split("__").str[0].map(_gram_stain)
                for stain in ("positive", "negative"):
                    sel = rows_t[gram == stain]
                    if len(sel):
                        put(f"{pre}_median_ba_tool_gram_{stain}",
                            float(sel["balanced_accuracy"].median()), "external_comparison",
                            "balanced_accuracy",
                            f"Median balanced accuracy of {tool}, Gram-{stain} organisms")
    lp = t.get("label_permutation_models")
    if lp is not None and "q" in lp:
        put("label_permutation.n_significant", int((lp["q"] < float(
            t.config["label_permutation"]["alpha"])).sum()), "label_permutation_models", "q",
            "Models whose label permutation q is below α")
        put("label_permutation.n_flagged", int((lp["q"] >= float(
            t.config["label_permutation"]["alpha"])).sum()), "label_permutation_models", "q",
            "Models whose q is at or above α: flagged, kept in every analysis (§8.6)")
    hyp = t.json("cross_model_dir", "hypotheses/hypotheses.json")
    if hyp is not None:
        _hypothesis_numbers(hyp, put, "hyp")
        put("hyp.n_organisms_without_reference", len(hyp.get("reference_missing") or []),
            "h7_determinants", "model_id",
            "Organisms without AMRFinderPlus calls, left out of H7 and H6")
        sens = hyp.get("sensitivity_label_permutation") or {}
        if sens.get("status") == "computed":
            put("lp_sensitivity.hyp.n_models_excluded", len(sens["excluded_models"]),
                "h1_models", "model_id", "Models the hypotheses are recomputed without (§14 item 9)")
            _hypothesis_numbers(sens, put, "lp_sensitivity.hyp", "Without the flagged models. ")
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

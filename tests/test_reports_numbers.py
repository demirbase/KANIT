#!/usr/bin/env python3
"""The aggregate keys of the numbers file (scripts/reports.py, numbers())."""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
pytestmark = pytest.mark.unit

DECISIONS = pd.DataFrame({
    "organism": ["ecoli", "ecoli", "kpneumoniae", "enterobacter_cloacae"],
    "antibiotic": ["ampicillin", "ciprofloxacin", "meropenem", "amikacin"],
    "drug_class": ["beta_lactams", "quinolones", "beta_lactams", "aminoglycosides"],
    "decision": ["included", "included", "included", "excluded"],
    "minority": [100, 200, 300, 50],
})
# Genomes that pass QC, per organism; the panel's own rule never sees these.
QC = {"ecoli": [True] * 10 + [False, False], "kpneumoniae": [True] * 5}
# Two models share genomes, so the union is smaller than the sum (3 + 3 + 2 = 8).
GENOMES = {("ecoli", "ampicillin"): ["g1", "g2", "g3"],
           ("ecoli", "ciprofloxacin"): ["g2", "g3", "g4"],
           ("kpneumoniae", "meropenem"): ["k1", "k2"]}
# The grade of the reported rule, per pattern: three of the six are weak or none.
PATTERNS = {("ecoli", "ampicillin"): ["weak", "confirmed"],
            ("ecoli", "ciprofloxacin"): ["none", "weak", "candidate"],
            ("kpneumoniae", "meropenem"): ["strong_novel"]}


# One win, one tie and one loss against the tool, spread over two organisms.
COMPARISON = {
    "ecoli": [("ecoli__ampicillin", "model", 1, 0.9),
              ("ecoli__ampicillin", "amrfinderplus", 1, 0.8),
              ("ecoli__ciprofloxacin", "model", 1, 0.8),
              ("ecoli__ciprofloxacin", "amrfinderplus", 1, 0.8)],
    "kpneumoniae": [("kpneumoniae__meropenem", "model", 1, 0.7),
                    ("kpneumoniae__meropenem", "amrfinderplus", 1, 0.9)],
}


class FakeTables:
    """Only the tables numbers() reads; everything else is absent."""

    config: dict = {"cv": {"n_bootstrap": 200}}

    def get(self, tid, organism=None, antibiotic=None):
        if tid == "panel_decisions":
            return DECISIONS
        key = (organism, antibiotic)
        if tid == "model_genomes" and key in GENOMES:
            ids = GENOMES[key]
            return pd.DataFrame({"Genome ID": ids, "label": [1] * len(ids)})
        if tid == "genome_qc" and organism in QC:
            return pd.DataFrame({"pass_overall": QC[organism]})
        if tid == "external_comparison" and organism in COMPARISON:
            return pd.DataFrame(COMPARISON[organism],
                                columns=["model_id", "tool", "assessable",
                                         "balanced_accuracy"])
        if tid == "grades_patterns" and key in PATTERNS:
            g = PATTERNS[key]
            return pd.DataFrame({"pattern_id": range(len(g)), "allele_aware_grade": g,
                                 "homolog_only_grade": g, "grade": g})
        return None


def _numbers():
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("amrtest_reports",
                                                  PROJECT_ROOT / "scripts" / "reports.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    df = m.numbers(FakeTables(), ["ecoli", "kpneumoniae", "enterobacter_cloacae"])
    return df.set_index("key")["value"].to_dict()


def test_panel_aggregates_count_only_included_pairs():
    """An organism whose every pair is excluded has no model, so it is not counted."""
    n = _numbers()
    assert n["panel.n_pairs"] == 4
    assert n["panel.n_included"] == 3
    assert n["panel.n_organisms_with_models"] == 2
    assert n["panel.n_antibiotics_included"] == 3
    assert n["panel.n_drug_classes_included"] == 2


def test_genomes_in_models_is_a_union_and_patterns_a_sum():
    """Genomes are shared between models; candidate patterns are model-specific."""
    n = _numbers()
    assert n["all.n_genomes_in_models"] == 6
    assert sum(n[f"{o}__{a}.n_genomes"] for o, a in GENOMES) == 8
    assert n["all.n_candidate_patterns"] == 6


def test_weak_or_none_share_is_measured_on_the_reported_rule():
    """Three of the six candidate patterns clear no evidence threshold."""
    n = _numbers()
    assert n["all.weak_or_none_share"] == 0.5


def test_the_share_is_absent_without_the_reported_grade():
    """An output carrying only the per-rule columns gives no share, not a zero."""
    class PerRuleOnly(FakeTables):
        def get(self, tid, organism=None, antibiotic=None):
            out = super().get(tid, organism, antibiotic)
            if tid == "grades_patterns" and out is not None:
                return out.drop(columns=["grade"])
            return out

    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("amrtest_reports3",
                                                  PROJECT_ROOT / "scripts" / "reports.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    df = m.numbers(PerRuleOnly(), ["ecoli", "kpneumoniae"])
    keys = set(df["key"])
    assert "all.weak_or_none_share" not in keys
    assert "all.n_candidate_patterns" in keys


def test_aggregates_are_absent_without_models():
    """Nothing is reported when no model has an output yet."""
    class Empty(FakeTables):
        def get(self, tid, organism=None, antibiotic=None):
            return DECISIONS if tid == "panel_decisions" else None

    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("amrtest_reports2",
                                                  PROJECT_ROOT / "scripts" / "reports.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    keys = set(m.numbers(Empty(), ["ecoli"])["key"])
    assert "all.n_genomes_in_models" not in keys
    assert "all.n_candidate_patterns" not in keys
    assert "all.weak_or_none_share" not in keys
    assert "panel.n_organisms_with_models" in keys


def test_the_sign_test_counts_pairs_but_tests_only_the_signed_ones():
    """The assessable count keeps the tie; the sign test drops it."""
    n = _numbers()
    assert n["compare.amrfinderplus_n_assessable"] == 3
    assert n["compare.amrfinderplus_n_model_better"] == 1
    assert n["compare.amrfinderplus_n_tool_better"] == 1
    assert n["compare.amrfinderplus_n_tied"] == 1
    assert n["compare.amrfinderplus_p"] == 1.0
    assert n["compare.amrfinderplus_median_ba_model"] == 0.8
    assert n["compare.amrfinderplus_median_ba_tool"] == 0.8
    # Two organisms carry the three pairs, so the interval is drawn over two clusters.
    assert n["compare.amrfinderplus_n_organisms"] == 2
    assert n["compare.amrfinderplus_median_diff"] == 0.0


def test_a_very_small_p_is_not_rounded_to_zero():
    """Four significant digits, so a p far below the fourth decimal survives."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("amrtest_reports4",
                                                  PROJECT_ROOT / "scripts" / "reports.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    assert m._four_significant(8.1861e-18) == 8.186e-18
    assert m._four_significant(0.058058) == 0.05806
    assert m._four_significant(0) == 0
    assert m._four_significant(None) is None


def test_grade_totals_sum_to_the_candidate_patterns():
    """The five grades partition the candidate set, so they add up to its size."""
    n = _numbers()
    grades = ("confirmed", "strong_novel", "candidate", "weak", "none")
    assert [n[f"all.n_{g}"] for g in grades] == [1, 1, 1, 2, 1]
    assert sum(n[f"all.n_{g}"] for g in grades) == n["all.n_candidate_patterns"]
    assert n["all.weak_or_none_share"] == 0.5


def test_per_organism_panel_counts_and_the_largest_minority():
    """An organism whose pairs all fall short keeps its counts and its minority."""
    n = _numbers()
    assert n["ecoli.n_pairs_assessed"] == 2
    assert n["ecoli.n_pairs_included"] == 2
    assert n["ecoli.max_minority"] == 200
    assert n["enterobacter_cloacae.n_pairs_assessed"] == 1
    assert n["enterobacter_cloacae.n_pairs_included"] == 0
    assert n["enterobacter_cloacae.max_minority"] == 50


def test_qc_genomes_are_summed_over_organisms():
    """Each organism has its own genomes, so the total is a sum, not a union."""
    n = _numbers()
    assert n["ecoli.n_genomes_qc_pass"] == 10
    assert n["all.n_genomes_qc_pass"] == 15

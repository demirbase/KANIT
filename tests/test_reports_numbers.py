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
})
# Two models share genomes, so the union is smaller than the sum (3 + 3 + 2 = 8).
GENOMES = {("ecoli", "ampicillin"): ["g1", "g2", "g3"],
           ("ecoli", "ciprofloxacin"): ["g2", "g3", "g4"],
           ("kpneumoniae", "meropenem"): ["k1", "k2"]}
PATTERNS = {("ecoli", "ampicillin"): 2, ("ecoli", "ciprofloxacin"): 3,
            ("kpneumoniae", "meropenem"): 1}


class FakeTables:
    """Only the tables numbers() reads; everything else is absent."""

    config: dict = {}

    def get(self, tid, organism=None, antibiotic=None):
        if tid == "panel_decisions":
            return DECISIONS
        key = (organism, antibiotic)
        if tid == "model_genomes" and key in GENOMES:
            ids = GENOMES[key]
            return pd.DataFrame({"Genome ID": ids, "label": [1] * len(ids)})
        if tid == "grades_patterns" and key in PATTERNS:
            n = PATTERNS[key]
            return pd.DataFrame({"pattern_id": range(n), "allele_aware_grade": ["weak"] * n,
                                 "homolog_only_grade": ["weak"] * n})
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
    assert "panel.n_organisms_with_models" in keys

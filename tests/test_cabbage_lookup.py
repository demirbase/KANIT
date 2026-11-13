"""KANIT lookup on the external isolates (lib/cabbage_lookup.py)."""
import math

import pandas as pd
import pytest
from lib import cabbage_lookup as cl


def _grades():
    return pd.DataFrame({
        "pattern_id": [1, 2, 3, 4, 5, 6],
        "grade": ["confirmed", "strong_novel", "candidate", "weak", "none", "candidate"],
        "allele_aware_card_state": ["b", "no_card_hit", "b", "no_card_hit", "no_card_hit",
                                    "card_hit_without_b"]})


def _prevalence():
    return pd.DataFrame({"pattern_id": [1, 2, 3, 4, 5, 6],
                         "direction": ["R", "R", "R", "S", "R", "S"]})


def test_markers_keep_graded_patterns_with_their_direction():
    mk = cl.markers(_grades(), _prevalence())
    assert sorted(mk["pattern_id"]) == [1, 2, 3, 4, 6]          # none is not a marker
    assert dict(zip(mk["pattern_id"], mk["direction"], strict=True))[4] == "S"
    with pytest.raises(ValueError):
        cl.markers(_grades(), _prevalence().iloc[:3])


def test_rules_take_the_r_markers_up_to_a_grade():
    mk = cl.markers(_grades(), _prevalence())
    assert cl.rule_patterns(mk, "upto_confirmed") == [1]
    assert cl.rule_patterns(mk, "upto_strong_novel") == [1, 2]
    assert cl.rule_patterns(mk, "upto_weak") == [1, 2, 3]      # 4 and 6 are S markers
    assert cl.rule_patterns(mk, "b1") == [1, 3]                 # catalogue part, any grade
    assert cl.rule_patterns(mk, "upto_weak", "S") == [4, 6]
    with pytest.raises(ValueError):
        cl.rule_patterns(mk, "upto_none")


def test_an_isolate_is_resistant_when_it_carries_a_marker_of_the_rule():
    mk = cl.markers(_grades(), _prevalence())
    presence = pd.DataFrame({"g1": [1, 0, 0, 0], "g2": [0, 1, 0, 1], "g3": [0, 0, 0, 1]},
                            index=[1, 2, 3, 4])                 # pattern 6 never called
    assert cl.predict(presence, mk, "upto_confirmed").tolist() == [1, 0, 0]
    assert cl.predict(presence, mk, "upto_strong_novel").tolist() == [1, 1, 0]
    assert cl.predict(presence, mk, "b1").tolist() == [1, 0, 0]
    assert cl.score(presence, mk).tolist() == [1.0, 0.0, -1.0]


def test_the_primary_rule_has_the_highest_pooled_balanced_accuracy():
    a = {"upto_confirmed": {"tp": 5, "fn": 5, "tn": 10, "fp": 0},
         "upto_strong_novel": {"tp": 8, "fn": 2, "tn": 9, "fp": 1},
         "upto_candidate": {"tp": 9, "fn": 1, "tn": 5, "fp": 5},
         "upto_weak": {"tp": 10, "fn": 0, "tn": 1, "fp": 9},
         "b1": {"tp": 10, "fn": 0, "tn": 10, "fp": 0}}           # b1 is never the primary
    pooled = cl.pooled([a, a])
    assert pooled["upto_confirmed"] == {"tp": 10, "fn": 10, "tn": 20, "fp": 0}
    assert cl.balanced_accuracy(pooled["upto_strong_novel"]) == pytest.approx(0.85)
    assert cl.select(pooled) == "upto_strong_novel"
    tie = {r: {"tp": 1, "fn": 1, "tn": 1, "fp": 1} for r in cl.UPTO}
    assert cl.select(tie) == "upto_confirmed"                   # a tie goes to the stricter
    assert math.isnan(cl.balanced_accuracy({"tp": 0, "fn": 0, "tn": 3, "fp": 1}))

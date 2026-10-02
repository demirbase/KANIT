#!/usr/bin/env python3
"""The panel rule (lib.panel): order of the rules, eligible genomes and counts."""
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import panel  # noqa: E402

pytestmark = pytest.mark.unit


def test_decide_rules_apply_in_order():
    # not a single drug: excluded whatever the counts
    assert panel.decide(900, 900, single_drug=False, registered=False,
                        min_minority=150) == (panel.EXCLUDED, "not a single drug")
    # below the threshold: excluded, registered or not
    assert panel.decide(149, 900, single_drug=True, registered=False,
                        min_minority=150) == (panel.EXCLUDED, "minority 149 < 150")
    # the threshold is inclusive, and the minority can be the susceptible class
    assert panel.decide(900, 150, single_drug=True, registered=True,
                        min_minority=150) == (panel.INCLUDED, "minority 150 >= 150")
    # passes the rule without a registry class: blocked, never silently included
    decision, _ = panel.decide(400, 400, single_drug=True, registered=False, min_minority=150)
    assert decision == panel.BLOCKED


def test_eligible_genomes_needs_qc_pass_and_a_lineage():
    qc = pd.DataFrame({"genome_id": ["g1", "g2", "g3", "g4"],
                       "pass_overall": [True, False, "True", None]})
    clusters = pd.DataFrame({"Genome ID": ["g1", "g2", "g3"], "Cluster": [1, 1, 2]})
    # g2 fails QC, g4 has a blank QC verdict, g3's "True" string still counts
    assert panel.eligible_genomes(qc, clusters) == {"g1", "g3"}
    # a genome that passes QC but has no lineage is not eligible
    assert panel.eligible_genomes(qc, clusters[clusters["Genome ID"] != "g3"]) == {"g1"}


def test_pair_rows_counts_only_eligible_tested_genomes():
    ph = pd.DataFrame({"Genome ID": ["g1", "g2", "g3", "g4"],
                       "ciprofloxacin": [1, 0, 1, None],
                       "extended spectrum beta lactamase": [1, 1, 0, 0]})
    rows = {r["antibiotic"]: r for r in panel.pair_rows("ecoli", ph, {"g1", "g2", "g4"},
                                                        min_minority=1)}
    cip = rows["ciprofloxacin"]
    assert (cip["n_tested"], cip["n_eligible"], cip["n_resistant"], cip["n_susceptible"]) == (3, 2, 1, 1)
    assert cip["decision"] == panel.INCLUDED and cip["drug_class"] == "quinolones"
    assert rows["extended spectrum beta lactamase"]["decision"] == panel.EXCLUDED


def test_pair_rows_rejects_values_other_than_0_and_1():
    ph = pd.DataFrame({"Genome ID": ["g1", "g2"], "ciprofloxacin": [1, 2]})
    with pytest.raises(ValueError, match="0 or 1"):
        panel.pair_rows("ecoli", ph, {"g1", "g2"}, min_minority=1)


def test_non_canonical_columns_are_reported():
    ph = pd.DataFrame({"Genome ID": ["g1"], "ciprofloxacin": [1], "fusidic acid": [0],
                       "ceftazidime/avibactam": [1], "mupirocin": [0]})
    # unknown names pass through unchanged and are left to the panel rules
    assert panel.non_canonical_columns(ph) == [("fusidic acid", "fusidic_acid"),
                                               ("ceftazidime/avibactam", "ceftazidime_avibactam")]


def test_included_pairs_reads_the_decisions(tmp_path):
    df = pd.DataFrame([
        {"organism": "ecoli", "antibiotic": "ciprofloxacin", "decision": panel.INCLUDED},
        {"organism": "ecoli", "antibiotic": "linezolid", "decision": panel.EXCLUDED},
    ])
    path = tmp_path / "panel_decisions.csv"
    df.to_csv(path, index=False)
    assert panel.included_pairs(path) == [("ecoli", "ciprofloxacin")]

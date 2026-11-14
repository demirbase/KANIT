#!/usr/bin/env python3
"""The keys the numbers file takes from the knowledge base (scripts/reports.py, _kb_numbers):
the mechanisms and gene families behind b = 1 and the NCBI context of the patterns."""
import importlib.util
import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from lib import knowledge_base as kb  # noqa: E402

pytestmark = pytest.mark.unit


def _reports():
    spec = importlib.util.spec_from_file_location("amrtest_reports_kb",
                                                  PROJECT_ROOT / "scripts" / "reports.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _df(**cols):
    return pd.DataFrame(cols)


@pytest.fixture
def db(tmp_path):
    """Two models. ecoli__ampicillin: pattern 1 confirmed (b = 1 on u1; u2 has a hit but no b),
    pattern 2 strong_novel towards S, pattern 3 from the association route (left out).
    enterococcus_faecium__vancomycin: pattern 1 confirmed, b = 1 on u4, which has no context."""
    path = tmp_path / "kanit.sqlite"
    con = sqlite3.connect(path)
    kb.create(con)
    con.execute("PRAGMA foreign_keys = OFF")
    kb.insert(con, "release", _df(kb_version=["1.0.0"], schema_version=["1.0"],
                                  created_at=["2026-11-14"], config_sha256=["x"], tools=["{}"]))
    models = ["ecoli__ampicillin", "enterococcus_faecium__vancomycin"]
    kb.insert(con, "model", _df(model_id=models, organism_id=["ecoli", "enterococcus_faecium"],
                                antibiotic_id=["ampicillin", "vancomycin"], n_genomes=[10, 10],
                                n_resistant=[5, 5], n_susceptible=[5, 5], n_lineages=[3, 3],
                                largest_lineage_share=[0.5, 0.5], n_unitigs=[9, 9],
                                n_patterns=[3, 1], evaluable=[1, 1]))
    m1, m2 = models
    kb.insert(con, "candidate", _df(model_id=[m1, m1, m1, m2], pattern_id=[1, 2, 3, 1],
                                    route=["model", "model", "association", "model"],
                                    source=["gain", "cpss", "pyseer", "gain"]))
    rows = [(m1, 1, "confirmed", "b"), (m1, 2, "strong_novel", "no_card_hit"),
            (m1, 3, "confirmed", "b"), (m2, 1, "confirmed", "b")]
    for mode, primary in (("allele_aware", 1), ("homolog_only", 0)):
        kb.insert(con, "grade", _df(
            model_id=[r[0] for r in rows], pattern_id=[r[1] for r in rows], mode=mode,
            primary_mode=primary, card_state=[r[3] for r in rows], card_reasons="",
            n_layers=4, layers_passed="p;m;c;y",
            # the secondary mode disagrees: it must not count
            grade=[r[2] if primary else "none" for r in rows]))
    kb.insert(con, "prevalence_result", _df(
        model_id=[m1, m1, m1, m2], pattern_id=[1, 2, 3, 1], present_resistant=4,
        present_susceptible=1, prev_resistant=0.8, prev_susceptible=0.2, delta=0.6,
        direction=["R", "S", "R", "R"], fisher_p=0.01, q=0.02, passes=1))
    kb.insert(con, "unitig", _df(unitig_id=["u1", "u2", "u3", "u4"],
                                 sequence=["A" * 40, "C" * 30, "G" * 50, "T" * 35],
                                 length=[40, 30, 50, 35]))
    kb.insert(con, "pattern_member", _df(model_id=[m1, m1, m1, m1, m2],
                                         pattern_id=[1, 1, 2, 3, 1],
                                         unitig_id=["u1", "u2", "u3", "u2", "u4"]))
    kb.insert(con, "card_annotation", _df(
        model_id=[m1, m1, m1, m2], unitig_id=["u1", "u2", "u3", "u4"], mode="allele_aware",
        state=["b", "card_hit_without_b", "no_card_hit", "b"], reasons="", n_b=1,
        located_genomes="g1"))
    kb.insert(con, "card_hit", _df(model_id=[m1, m1, m2, m2],
                                   unitig_id=["u1", "u2", "u4", "u4"],
                                   aro_accession=["A1", "A2", "A3", "A1"]))
    kb.insert(con, "aro", _df(
        aro_accession=["A1", "A2", "A3"], name=["TEM-1", "pump", "vanA"], model_type="homolog",
        gene_family=["TEM beta-lactamase", "MFS pump",
                     "glycopeptide resistance gene cluster; van ligase"],
        mechanism=["antibiotic inactivation", "antibiotic efflux",
                   "antibiotic target alteration"]))
    kb.insert(con, "unitig_context", _df(
        unitig_id=["u1", "u2", "u3"], organism_id="ecoli", source="ncbi_nt_remote",
        queried_on="2026-11-12", n_hits=[50, 50, 10], plasmid_share=[0.8, 0.2, 0.0],
        gene=["blaTEM", "x", None]))
    con.commit()
    con.close()
    return path


def _keys(db):
    out = {}

    def put(key, value, table, column, description, digits="auto"):
        out[key] = value
    _reports()._kb_numbers(db, put)
    return out


def test_release_and_table_counts(db):
    k = _keys(db)
    assert k["kb.version"] == "1.0.0"
    assert k["kb.n_rows.candidate"] == 4
    assert k["kb.n_views"] == 4


def test_mechanisms_count_only_the_aros_on_b_unitigs(db):
    k = _keys(db)
    assert k["kb.card.n_patterns_b"] == 2           # the association pattern is left out
    assert k["kb.card.n_aros"] == 2                 # A2 sits on u2, which has no b
    assert k["kb.card.n_mechanisms"] == 2
    assert k["kb.card.n_gene_families"] == 3        # A3's two families count apart
    assert k["kb.card.n_gene_families_multi_organism"] == 1      # TEM, in both
    assert k["kb.card.ecoli.share_inactivation"] == 1.0
    assert "kb.card.ecoli.share_efflux" not in k
    assert k["kb.card.enterococcus_faecium.share_target_alteration"] == 1.0
    assert k["kb.card.enterococcus_faecium.share_inactivation"] == 1.0


def test_context_takes_the_longest_searched_member(db):
    k = _keys(db)
    assert k["kb.context.n_patterns"] == 3
    assert k["kb.context.n_patterns_searched"] == 2          # u4 was not searched
    assert k["kb.context.confirmed.n_patterns"] == 2
    assert k["kb.context.confirmed.n_with_hit"] == 1
    assert k["kb.context.confirmed.median_plasmid_share"] == 0.8      # u1, not u2
    assert k["kb.context.confirmed.share_mostly_plasmid"] == 1.0
    assert k["kb.context.strong_novel.S.n_patterns"] == 1
    assert k["kb.context.strong_novel.S.share_no_plasmid"] == 1.0
    assert k["kb.context.strong_novel.S.share_with_gene"] == 0.0
    assert k["kb.context.strong_novel.R.n_patterns"] == 0
    assert "kb.context.strong_novel.R.median_plasmid_share" not in k


# ---- the three case studies (_kb_case_numbers) ----

def _pattern(con, model, organism, pid, grade, direction, pi, gain, delta, unitigs, carriers,
             pyseer_p=1e-5, mda=0.001, hits=50, card_state="no_card_hit"):
    """One graded pattern of the model route with all its layers. unitigs: {id: sequence};
    carriers: one character per genome, '1' for a carrier."""
    key = dict(model_id=[model], pattern_id=[pid])
    kb.insert(con, "candidate", _df(**key, route="model", source="cpss" if gain is None else "gain",
                                    gain_rank=None, total_gain=[gain]))
    for mode, primary in (("allele_aware", 1), ("homolog_only", 0)):
        kb.insert(con, "grade", _df(**key, mode=mode, primary_mode=primary, card_state=card_state,
                                    card_reasons="", n_layers=4, layers_passed="p;m;c;y",
                                    grade=grade))
    n = len(carriers)
    kb.insert(con, "prevalence_result", _df(**key, present_resistant=3, present_susceptible=1,
                                            prev_resistant=0.75, prev_susceptible=0.25,
                                            delta=delta, direction=direction, fisher_p=0.01,
                                            q=0.02, passes=1))
    kb.insert(con, "mda_result", _df(**key, n_fold_models_using=5, auc_observed=0.9,
                                     auc_permuted=0.8, mda=mda, n_permuted_ge_observed=0,
                                     p=0.001, q=0.01, passes=1))
    kb.insert(con, "cpss_result", _df(**key, in_prefilter=1, chi2=10.0, n_selected=200, pi=pi,
                                      passes=1))
    kb.insert(con, "pyseer_result", _df(**key, af=0.5, beta=1.0, beta_se=0.1, lrt_p=pyseer_p,
                                        notes="", passes=1))
    blob = int(carriers.ljust(-(-n // 8) * 8, "0"), 2).to_bytes(-(-n // 8), "big")
    kb.insert(con, "pattern", _df(**key, n_members=len(unitigs), n_present=carriers.count("1"),
                                  carriers=[blob]))
    for uid, seq in unitigs.items():
        kb.insert(con, "unitig", _df(unitig_id=[uid], sequence=[seq], length=[len(seq)]))
        kb.insert(con, "pattern_member", _df(**key | {"unitig_id": [uid]}))
        kb.insert(con, "card_annotation", _df(model_id=[model], unitig_id=[uid],
                                              mode="allele_aware", state=card_state, reasons="",
                                              n_b=0, located_genomes=""))
    longest = max(unitigs, key=lambda u: len(unitigs[u]))
    if hits is not None:
        kb.insert(con, "unitig_context", _df(unitig_id=[longest], organism_id=[organism],
                                             source="ncbi_nt_remote", queried_on="2026-11-12",
                                             n_hits=[hits], plasmid_share=[0.0], gene=[None]))


@pytest.fixture
def cases_db(tmp_path):
    path = tmp_path / "kanit.sqlite"
    con = sqlite3.connect(path)
    kb.create(con)
    con.execute("PRAGMA foreign_keys = OFF")
    models = {"ecoli__amp": ("ecoli", 1), "kpneumoniae__amp": ("kpneumoniae", 0),
              "efaecium__tet": ("efaecium", 1), "saureus__cip": ("saureus", 1)}
    for m, (org, amrfinder) in models.items():
        kb.insert(con, "model", _df(model_id=[m], organism_id=[org], antibiotic_id=["x"],
                                    n_genomes=[8], n_resistant=[4], n_susceptible=[4],
                                    n_lineages=[3], largest_lineage_share=[0.5], n_unitigs=[9],
                                    n_patterns=[3], evaluable=[1]))
        kb.insert(con, "external_comparison", _df(model_id=[m], tool=["amrfinderplus"],
                                                  assessable=[amrfinder]))
        kb.insert(con, "model_metric", _df(model_id=[m], arm=["lineage_aware"],
                                           metric=["roc_auc"], value=[0.9]))
    s = "ACGTTGCAAC"
    # case 1: the larger gain wins; kpneumoniae is not assessed by AMRFinderPlus
    _pattern(con, "ecoli__amp", "ecoli", 1, "confirmed", "R", 1.0, 10.0, 0.5, {"E1": "A" * 31},
             "11110000")
    _pattern(con, "ecoli__amp", "ecoli", 2, "confirmed", "R", 1.0, 50.0, 0.5, {"E2": "C" * 31},
             "11100000")
    _pattern(con, "kpneumoniae__amp", "kpneumoniae", 1, "confirmed", "R", 1.0, 999.0, 0.5,
             {"K1": "G" * 31}, "11110000")
    # case 2: ecoli is taken; π and gain tie, |Δ| decides; a pattern without a hit is left out
    _pattern(con, "ecoli__amp", "ecoli", 3, "strong_novel", "R", 1.0, 100.0, 0.9,
             {"E3": "AC" * 16}, "11110000")
    _pattern(con, "kpneumoniae__amp", "kpneumoniae", 2, "strong_novel", "R", 1.0, None, 0.3,
             {"K2": "AG" * 16}, "11000000", pyseer_p=1e-3, mda=0.0001)
    _pattern(con, "efaecium__tet", "efaecium", 1, "strong_novel", "R", 1.0, None, 0.9,
             {"F1": "AT" * 16}, "11110000", pyseer_p=1e-9, mda=0.01)
    _pattern(con, "efaecium__tet", "efaecium", 2, "strong_novel", "R", 1.0, None, 0.95,
             {"F2": "CG" * 16}, "11110000", hits=0)
    _pattern(con, "efaecium__tet", "efaecium", 3, "confirmed", "R", 1.0, 30.0, 0.8,
             {"F3": "GT" * 16}, "11111000")
    # case 3: towards S, from neither ecoli nor efaecium; one base from a confirmed unitig
    _pattern(con, "saureus__cip", "saureus", 1, "confirmed", "R", 1.0, 40.0, 0.8, {"S1": s},
             "11110000")
    _pattern(con, "saureus__cip", "saureus", 2, "strong_novel", "S", 1.0, 20.0, -0.8,
             {"S2": "ACGTCGCAAC"[::-1].translate(str.maketrans("ACGT", "TGCA"))}, "00001111")
    _pattern(con, "kpneumoniae__amp", "kpneumoniae", 3, "strong_novel", "S", 0.9, 500.0, -0.9,
             {"K3": "TC" * 16}, "00001111")
    con.commit()
    con.close()
    return path


def test_cases_follow_the_rule(cases_db):
    k = _keys(cases_db)
    assert k["kb.case1.model_id"] == "ecoli__amp" and k["kb.case1.pattern_id"] == "2"
    assert k["kb.case1.n_eligible"] == 4 and k["kb.case1.n_tied"] == 1
    assert k["kb.case1.n_carriers"] == 3 and k["kb.case1.n_carriers_with_confirmed"] == 3
    assert k["kb.case2.model_id"] == "efaecium__tet" and k["kb.case2.pattern_id"] == "1"
    assert k["kb.case2.n_eligible"] == 2                # ecoli taken, F2 has no hit
    assert k["kb.case2.n_tied"] == 2 and k["kb.case2.alternatives_agree"] == 1
    assert k["kb.case2.n_carriers_with_confirmed"] == 4
    assert k["kb.case3.model_id"] == "saureus__cip" and k["kb.case3.unitig_id"] == "S2"
    assert k["kb.case3.neighbour_pattern_id"] == "1"   # S2 is S1 with one base changed, other strand
    assert k["kb.case3.n_carriers_with_neighbour"] == 0
    assert k["kb.case3.n_carriers"] == 4 and k["kb.case3.neighbour_n_carriers"] == 4
    assert "kb.case1.neighbour_pattern_id" not in k

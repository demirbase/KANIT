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

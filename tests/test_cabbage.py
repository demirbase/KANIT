#!/usr/bin/env python3
"""External validation on CABBAGE: the eligible isolates (scripts/lib/cabbage.py)."""
import sys
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import cabbage, registry  # noqa: E402

pytestmark = pytest.mark.unit


def _row(bs, species, ab, std="CLSI", pheno="resistant", ecoff="no", asm=None, sra=None,
         db="NCBI_antibiogram"):
    return {"BioSample_ID": bs, "SRA_accession": sra, "assembly_ID": asm, "species": species,
            "database": db, "antibiotic_name": ab, "ast_standard": std,
            "resistance_phenotype": pheno, "used_ECOFF": ecoff}


EC, KP = "Escherichia coli", "Klebsiella pneumoniae"
ROWS = [
    _row("S1", EC, "ciprofloxacin", asm="GCA_000000001.1"),
    _row("S1", EC, "ciprofloxacin", std="EUCAST", asm="GCA_000000001.1"),
    _row("S2", EC, "ciprofloxacin", pheno="susceptible", asm="GCA_000000002.1"),   # tie
    _row("S2", EC, "ciprofloxacin", asm="GCA_000000002.1"),
    _row("S2", EC, "ampicillin", pheno="susceptible", asm="GCA_000000002.1"),
    _row("S3", "Enterobacter hormaechei", "trimethoprim-sulfamethoxazole", pheno="susceptible",
         asm="GCA_000000003.1"),
    _row("S3", "Enterobacter hormaechei", "trimethoprim-sulfamethoxazole", pheno="susceptible",
         asm="GCA_000000003.2"),
    _row("S4", "Enterobacter kobei", "ciprofloxacin", asm="GCA_000000004.1"),     # species
    _row("S5", EC, "ciprofloxacin", std="NCCLS", asm="GCA_000000005.1"),         # standard
    _row("S6", EC, "ciprofloxacin", pheno="intermediate", asm="GCA_000000006.1"),
    _row("S7", EC, "ciprofloxacin", ecoff="yes", asm="GCA_000000007.1"),
    _row("S8", EC, "ciprofloxacin", asm="GCA_000000008.1"),
    _row("S8", EC, "ampicillin", std="NARMS", db="PATRIC;NARMS", asm="GCA_000000008.1"),
    _row("S9", EC, "ciprofloxacin", asm="GCA_000000009.1"),                      # BioSample
    _row("S10", EC, "ciprofloxacin", asm="GCA_000000010.1"),                     # assembly
    _row("S11", EC, "ciprofloxacin"),                                            # no assembly
    _row("S12", EC, "ciprofloxacin", asm="GCA_000000012.1"),
    _row("S12", EC, "ciprofloxacin", asm="GCA_000000013.1"),                     # two assemblies
    _row("S13", KP, "amoxicillin-clavulanic acid", asm="GCA_000000014.1"),
    _row("S13", KP, "cefpodoxime-clavulanic acid", pheno="susceptible", asm="GCA_000000014.1"),
]
IDS = {"biosample": {"S9"}, "sra": set(), "assembly": {"GCA_000000010"}}


def test_names_and_accessions():
    assert cabbage.antibiotic("trimethoprim-sulfamethoxazole") == "trimethoprim_sulfamethoxazole"
    assert cabbage.antibiotic("cefpodoxime-clavulanic acid") == "cefpodoxime_clavulanic_acid"
    assert cabbage.antibiotic("Imipenem") == "imipenem" and cabbage.antibiotic(" ") is None
    assert cabbage.accessions(["SRR1, SRR2", None, "srr3"]) == {"SRR1", "SRR2", "SRR3"}
    assert cabbage.unversioned("GCA_000001405.15") == "GCA_000001405"
    g = pd.DataFrame({"biosample_accession": ["SAMN1", ""], "sra_accession": ["SRR1,SRR2", ""],
                      "assembly_accession": ["GCA_000000001.2", ""]})
    assert cabbage.bvbrc_ids(g) == {"biosample": {"SAMN1"}, "sra": {"SRR1", "SRR2"},
                                    "assembly": {"GCA_000000001"}}
    idx = cabbage.species_index(registry.load_organisms())
    assert idx["Enterobacter bugandensis"] == "enterobacter_cloacae" and "Enterobacter kobei" not in idx


def test_select_follows_every_rule():
    isolates, phenotypes, counts = cabbage.select(pd.DataFrame(ROWS), registry.load_organisms(), IDS)
    assert isolates["biosample_id"].tolist() == ["S1", "S2", "S3", "S13"]
    assert dict(zip(isolates["biosample_id"], isolates["assembly_id"], strict=True))["S3"] == \
        "GCA_000000003.2"                                      # the latest version
    p = {(r.biosample_id, r.antibiotic): r for r in phenotypes.itertuples()}
    assert p[("S1", "ciprofloxacin")].label == 1 and p[("S1", "ciprofloxacin")].n_records == 2
    assert ("S2", "ciprofloxacin") not in p and p[("S2", "ampicillin")].label == 0
    assert p[("S3", "trimethoprim_sulfamethoxazole")].label == 0
    assert {a for b, a in p if b == "S13"} == {"amoxicillin_clavulanic_acid",
                                               "cefpodoxime_clavulanic_acid"}
    ec = counts["ecoli"]
    assert [ec[s] for s in cabbage.STEPS] == [10, 9, 8, 7, 7, 6, 4, 2]
    assert ec["cells_tied"] == 1 and ec["cells"] == 2
    assert counts["enterobacter_cloacae"]["assembly"] == 1 and counts["kpneumoniae"]["cells"] == 2
    t = cabbage.pair_counts(phenotypes).set_index(["organism", "antibiotic"])
    assert t.loc[("ecoli", "ciprofloxacin"), "n_resistant"] == 1
    assert t.loc[("ecoli", "ampicillin"), "n_susceptible"] == 1


def test_select_needs_cabbage_columns():
    with pytest.raises(ValueError, match="columns missing"):
        cabbage.select(pd.DataFrame(ROWS).drop(columns=["used_ECOFF"]),
                       registry.load_organisms(), IDS)

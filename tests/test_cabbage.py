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


# ---- NCBI Datasets client and the fetch step ------------------------------------------

def _zip(acc: str, fasta: bytes) -> bytes:
    import io
    import zipfile
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr(f"ncbi_dataset/data/{acc}/{acc}_ASM1_genomic.fna", fasta)
        z.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}")
    return buf.getvalue()


class FakeNcbi:
    """fetch(url) of the Datasets client: genome reports, lineages and packages."""

    def __init__(self, genomes, taxa):
        self.genomes, self.taxa, self.calls = genomes, taxa, []

    def __call__(self, url: str) -> bytes:
        import json
        import urllib.parse
        self.calls.append(url)
        path = urllib.parse.unquote(url.split("/v2/", 1)[1])
        if path.startswith("genome/accession/") and "/dataset_report" in path:
            accs = path.split("/")[2].split(",")
            return json.dumps({"reports": [
                {"accession": a, "organism": {"tax_id": g["tax"], "organism_name": g["name"]},
                 "assembly_stats": {"number_of_contigs": g["contigs"],
                                    "total_sequence_length": str(g["length"])}}
                for a in accs if (g := self.genomes.get(a))]}).encode()
        if path.startswith("taxonomy/taxon/"):
            ids = [int(t) for t in path.split("/")[2].split(",")]
            return json.dumps({"reports": [{"taxonomy": {"tax_id": t, "parents": self.taxa[t]}}
                                           for t in ids if t in self.taxa]}).encode()
        acc = path.split("/")[2]
        return _zip(acc, self.genomes[acc]["fasta"])


def test_datasets_client(monkeypatch):
    import urllib.error

    from lib import ncbi_datasets as nd
    fake = FakeNcbi({"GCA_1.1": {"tax": 562, "name": "E. coli", "contigs": 2, "length": 8,
                                 "fasta": b">a\nACGT\n>b\nACGT\n"},
                     "GCA_2.1": {"tax": 562, "name": "E. coli", "contigs": 1, "length": 4,
                                 "fasta": b">a\nACGT\n"}}, {562: [1, 561]})
    monkeypatch.setattr(nd, "BATCH", 1)
    slept = []
    api = nd.Datasets(fake, sleep=slept.append)
    assert set(api.genome_reports(["GCA_1.1", "GCA_2.1", "GCA_9.1"])) == {"GCA_1.1", "GCA_2.1"}
    assert sum("dataset_report" in c for c in fake.calls) == 3          # one per batch
    assert api.lineages([562]) == {562: {562, 561, 1}}
    assert api.fasta("GCA_1.1") == b">a\nACGT\n>b\nACGT\n"
    rep = api.genome_reports(["GCA_1.1"])["GCA_1.1"]
    assert nd.report_numbers(rep) == ("2", "8", 562, "E. coli")
    tries = iter([urllib.error.HTTPError("u", 429, "busy", {}, None), b'{"reports": []}'])

    def flaky(url):
        r = next(tries)
        if isinstance(r, Exception):
            raise r
        return r
    assert nd.Datasets(flaky, sleep=slept.append, wait=1).lineages([1]) == {}
    assert 1 in slept                                                   # waited once


def test_fetch_checks_taxon_and_assembly(tmp_path):
    import importlib.util

    from lib import ncbi_datasets as nd
    from lib.config import load_config
    spec = importlib.util.spec_from_file_location("amrtest_19", PROJECT_ROOT / "scripts" / "19_cabbage.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    config = load_config()
    config["paths_organism"] = {**config["paths_organism"], "cabbage_dir": str(tmp_path / "out"),
                                "cabbage_genomes_dir": str(tmp_path / "g" / "{organism}")}
    (tmp_path / "out").mkdir()
    pd.DataFrame({"organism": ["ecoli"] * 4, "biosample_id": ["S1", "S2", "S3", "S4"],
                  "assembly_id": ["GCA_1.1", "GCA_2.1", "GCA_3.1", "GCA_4.1"],
                  "sra_accession": "", "sources": "NCBI_antibiogram"}).to_csv(
        tmp_path / "out" / "cabbage_isolates.csv", index=False)
    two = b">a\nACGT\n>b\nACGT\n"
    fake = FakeNcbi({"GCA_1.1": {"tax": 562, "name": "Escherichia coli", "contigs": 2, "length": 8, "fasta": two},
                     "GCA_2.1": {"tax": 1313, "name": "Streptococcus pneumoniae", "contigs": 2, "length": 8, "fasta": two},
                     "GCA_3.1": {"tax": 562, "name": "Escherichia coli", "contigs": 3, "length": 8, "fasta": two}},
                    {562: [1, 561], 1313: [1, 1301]})
    api = nd.Datasets(fake, sleep=lambda s: None)
    r = m.fetch(config, workers=2, api=api).set_index("biosample_id")
    assert r.loc["S1", "status"] == "passed" and r.loc["S1", "bytes"] == len(two)
    assert r.loc["S2", "status"] == "excluded" and "1313" in r.loc["S2", "problem"]
    assert r.loc["S3", "status"] == "failed" and "the record has 3" in r.loc["S3", "problem"]
    assert r.loc["S4", "status"] == "failed" and r.loc["S4", "problem"] == "no NCBI genome report"
    assert (tmp_path / "g" / "ecoli" / "S1.fna").read_bytes() == two
    n_downloads = sum("/download" in c for c in fake.calls)
    m.fetch(config, workers=2, api=api)                                # S1 is not fetched again
    assert sum("/download" in c for c in fake.calls) == n_downloads + 1   # only S3 again

#!/usr/bin/env python3
"""BV-BRC snapshot (lib.bvbrc, 00a_download_bvbrc.py): name normalisation, the
phenotype rules of protocol §2.2 and the snapshot against a stand-in API."""
import importlib.util
import json
import re
import sys
import urllib.error
import urllib.parse
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import bvbrc, registry  # noqa: E402
from lib.bvbrc import clean_amr_table, pivot_binary, resolve  # noqa: E402

pytestmark = pytest.mark.unit


# ---- antibiotic name normalisation (registry) -----------------------------------------
def test_normalize_aliases_and_typos():
    n = registry.normalize_antibiotic
    assert n("amipicillin_sulbactam") == "ampicillin_sulbactam"
    assert n("ampicillin/sulbactam") == "ampicillin_sulbactam"     # BV-BRC spelling -> canonical
    assert n("rifampicin") == "rifampin"
    assert n("cefalotin") == "cephalothin"
    assert n("tigecyklin") == "tigecycline"
    assert n("amoxicillin_clavulanat") == "amoxicillin_clavulanic_acid"
    assert n("co-trimoxazole") == "trimethoprim_sulfamethoxazole"
    assert n("sulfamethoxazole/trimethoprim") == "trimethoprim_sulfamethoxazole"
    assert n("  GENTAMICIN ") == "gentamicin"
    assert n("some_unlisted_drug") == "some_unlisted_drug"   # unknown kept, not dropped
    assert n(None) is None


# ---- phenotype rules (§2.2) -----------------------------------------------------------
LAB = "Laboratory Method"


def _records(rows):
    return pd.DataFrame(rows, columns=["genome_id", "antibiotic", "resistant_phenotype",
                                       "evidence", "testing_standard", "testing_standard_year"])


def test_record_filters():
    df = _records([
        ("1", "ampicillin", "Resistant", LAB, "EUCAST", "2020"),
        ("1", "ampicillin", "Resistant", LAB, "EUCAST", "2020"),         # identical: counts once
        ("1", "tetracycline", "Resistant", LAB, "NARMS", "2021"),        # other standard
        ("1", "cefoxitin", "Resistant", LAB, "", ""),                    # no standard
        ("1", "ciprofloxacin", "Intermediate", LAB, "CLSI", "2021"),     # not R/S
        ("1", "gentamicin", "Susceptible", LAB, "EUCAST, CLSI", "2019"),  # combined form kept
        ("1", "meropenem", "Susceptible", LAB, "clsi", "2019"),          # case-insensitive
        ("2", "meropenem", "Resistant", "", "CLSI", "2019"),             # evidence empty
        ("2", "imipenem", "Resistant", "Computational Method", "CLSI", "2019"),
        ("2", "fluoroquinolones", "Resistant", LAB, "CLSI", "2019"),     # not a registered drug
    ])
    cleaned, rep = clean_amr_table(df)
    assert {(r.genome_id, r.antibiotic): r.label for r in cleaned.itertuples()} == {
        ("1", "ampicillin"): 1, ("1", "gentamicin"): 0, ("1", "meropenem"): 0,
        ("2", "fluoroquinolones"): 1}
    assert (rep["records"], rep["records_distinct"], rep["records_laboratory"],
            rep["records_eucast_clsi"], rep["records_resistant_susceptible"]) == (10, 9, 7, 5, 4)
    assert rep["phenotypes_dropped"] == {"intermediate": 1}
    assert rep["antibiotics_not_registered"] == ["fluoroquinolones"]


def test_conflicts():
    assert resolve([1, 1, 0], [2019, 2020, 2018]) == (1, "majority")
    assert resolve([1, 0], [2015, 2022]) == (0, "recent_year")
    assert resolve([1, 0], [2015, np.nan]) == (1, "recent_year")     # the year that exists
    assert resolve([1, 0], [np.nan, np.nan]) == (None, "tied")
    assert resolve([1, 0, 1, 0], [2020, 2020, 2019, 2019]) == (None, "tied")  # tied in 2020
    assert resolve([1, 0, 0, 1], [2021, 2020, 2020, 2019]) == (1, "recent_year")
    assert resolve([0, 0], ["", ""]) == (0, "single")
    df = _records([("3", "gentamicin", "Resistant", LAB, "EUCAST", "2015"),
                   ("3", "gentamicin", "Susceptible", LAB, "CLSI", "2015"),
                   ("3", "ampicillin", "Resistant", LAB, "EUCAST", "2015"),
                   ("3", "ampicillin", "Resistant", LAB, "CLSI", "2016"),
                   ("3", "ampicillin", "Susceptible", LAB, "CLSI", "2017")])
    cleaned, rep = clean_amr_table(df)
    assert cleaned.to_dict("records") == [{"genome_id": "3", "antibiotic": "ampicillin",
                                           "label": 1}]
    assert (rep["cells"], rep["cells_conflicting"], rep["cells_resolved_by_majority"],
            rep["cells_dropped_tied"]) == (2, 2, 1, 1)


def test_pivot_binary():
    wide = pivot_binary(pd.DataFrame({"genome_id": ["1", "1", "2"],
                                      "antibiotic": ["ampicillin", "gentamicin", "ampicillin"],
                                      "label": [1, 0, 1]}))
    assert list(wide.columns) == ["Genome ID", "ampicillin", "gentamicin"]
    assert pd.isna(wide.loc[wide["Genome ID"] == "2", "gentamicin"]).all()


def test_check_fasta():
    fa = b">c1 x\nACGT\nAC\n>c2\nGGG\n"
    assert bvbrc.check_fasta(fa, "2", "9") == (2, 9, "")
    assert bvbrc.check_fasta(fa, "3", "9")[2] == "2 sequences, the record has 3"
    assert bvbrc.check_fasta(fa, "", "") == (2, 9, "")
    assert bvbrc.check_fasta(b"<html>", "1", "4")[2] == "not FASTA"


# ---- a stand-in for the BV-BRC API ------------------------------------------------------
GENOMES = pd.DataFrame({
    "genome_id": ["562.1", "562.2", "562.3", "562.4", "562.5"],
    "genome_name": ["E. coli a", "E. coli b", "E. coli c", "E. coli d", "E. coli e"],
    "taxon_id": ["562", "83334", "562", "562", "562"],
    "assembly_accession": ["GCA_1.1", "", "GCA_3.1", "", ""],
    "sra_accession": ["SRR1", "SRR2,SRR3", "", "SRR4", ""],
    "biosample_accession": ["SAMN1", "SAMN2", "SAMN3", "SAMN4", ""],
    "contigs": ["2", "1", "1", "1", "1"], "genome_length": ["9", "4", "4", "4", "4"]})
RECORDS = pd.DataFrame(
    [("562.1", "ampicillin", "Resistant", LAB, "CLSI", "2020"),
     ("562.1", "ciprofloxacin", "Susceptible", LAB, "EUCAST", "2021"),
     ("562.2", "ampicillin", "Susceptible", LAB, "CLSI", "2019"),
     ("562.3", "ampicillin", "Resistant", LAB, "CLSI", "2019"),
     ("562.4", "ampicillin", "Resistant", LAB, "EUCAST", "2018")],
    columns=["genome_id", "antibiotic", "resistant_phenotype", "evidence", "testing_standard",
             "testing_standard_year"])                            # 562.5: no records
FASTA = {"562.1": b">a\nACGT\n>b\nACGTA\n", "562.2": b">a\nACGT\n", "562.3": b">a\nTTTT\n",
         "562.4": b">a\nACG\n"}                                   # 562.4: shorter than its record


class FakeApi:
    def __init__(self, fail_first: int = 0):
        self.calls, self.fail_first = [], fail_first

    def __call__(self, url, data, accept):
        self.calls.append(url)
        if self.fail_first:
            self.fail_first -= 1
            raise urllib.error.HTTPError(url, 503, "busy", {}, None)
        if data is None and url.endswith("/api/"):
            return b"<div>Version: 1.9.3</div>", {}
        if data is None:
            gid = re.search(r"eq\(genome_id,([^)]+)\)", url).group(1)
            return FASTA[gid], {}
        rql = urllib.parse.unquote(data.decode())
        if url.endswith("/genome/"):
            assert "in(taxon_lineage_ids,(562))" in rql and "eq(public,true)" in rql
            table = GENOMES
        else:
            ids = re.search(r"in\(genome_id,\(([^)]*)\)\)", rql).group(1).split(",")
            assert f'eq(evidence,"{LAB}")' in rql
            table = RECORDS[RECORDS["genome_id"].isin(ids)]
        fields = re.search(r"select\(([^)]*)\)", rql).group(1).split(",")
        size, start = map(int, re.search(r"limit\((\d+),(\d+)\)", rql).groups())
        part = table.reindex(columns=fields, fill_value="").iloc[start:start + size]
        return (part.to_csv(index=False).encode(),
                {"Content-Range": f"items {start}-{start + len(part)}/{len(table)}"})


def test_api_retries_and_paging(monkeypatch):
    monkeypatch.setattr(bvbrc, "PAGE", 2)
    fake = FakeApi(fail_first=1)
    api = bvbrc.Api(fake, sleep=lambda s: None)
    g = bvbrc.genomes(api, [562])
    assert g["genome_id"].tolist() == GENOMES["genome_id"].tolist() and len(fake.calls) == 4
    assert api.version() == "1.9.3"

    def bad(url, data, accept):
        raise urllib.error.HTTPError(url, 400, "bad query", {}, None)
    with pytest.raises(urllib.error.HTTPError):
        bvbrc.Api(bad, sleep=lambda s: None).version()


def _script():
    spec = importlib.util.spec_from_file_location("amrtest_00a", PROJECT_ROOT / "scripts" /
                                                  "00a_download_bvbrc.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_snapshot_end_to_end(tmp_path, monkeypatch):
    monkeypatch.setattr(bvbrc, "PAGE", 2)
    monkeypatch.setattr(bvbrc, "ID_BATCH", 2)
    m = _script()
    meta, genomes_dir = tmp_path / "meta", tmp_path / "genomes"
    api = bvbrc.Api(FakeApi(), sleep=lambda s: None)
    q = m.query("ecoli", meta, api)
    assert (q["api_version"], q["n_genomes_at_or_below"], q["n_genomes_with_records"],
            q["n_genomes_phenotyped"]) == ("1.9.3", 5, 4, 4)
    g = pd.read_csv(meta / "genomes.csv", dtype=str, keep_default_na=False)
    assert g.loc[g["genome_id"] == "562.2", "sra_accession"].item() == "SRR2,SRR3"
    report = m.fetch(meta, genomes_dir, api, workers=2)
    assert report.set_index("genome_id")["status"].to_dict() == {
        "562.1": "passed", "562.2": "passed", "562.3": "passed", "562.4": "failed"}
    assert not (genomes_dir / "562.4.fna").exists()
    fake = FakeApi()
    again = m.fetch(meta, genomes_dir, bvbrc.Api(fake, sleep=lambda s: None), workers=1)
    assert len(fake.calls) == 1 and again["status"].tolist() == report["status"].tolist()
    snap = m.freeze(meta, meta / "amr_phenotypes.csv")
    assert (snap["n_genomes"], snap["n_assemblies_failed"], snap["n_antibiotics"]) == (3, 1, 2)
    wide = pd.read_csv(meta / "amr_phenotypes.csv", dtype={"Genome ID": str})
    assert wide["Genome ID"].tolist() == ["562.1", "562.2", "562.3"]
    assert json.loads((meta / "snapshot.json").read_text())["query"]["taxids"] == [562]
    assert m.verify(meta, genomes_dir) == []
    (genomes_dir / "562.2.fna").write_bytes(b">a\nAAAA\n")
    assert [Path(p).name for p in m.verify(meta, genomes_dir)] == ["562.2.fna"]
    with pytest.raises(SystemExit):
        m.query("ecoli", meta, api)                         # frozen: not queried again

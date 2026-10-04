#!/usr/bin/env python3
"""NCBI context (lib.ncbi_context, 18_ncbi_context.py) with a fake NCBI: no network."""
import importlib.util
import io
import json
import sys
import urllib.parse
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import knowledge_base as kb  # noqa: E402
from lib import ncbi_context as nc  # noqa: E402

pytestmark = pytest.mark.unit

TEXT = "BLASTN 2.17\n...\n  Database: Nucleotide collection (nt)\n    Posted date:  Oct 1, 2026\n"


def _xml(queries) -> str:
    its = []
    for q, hits in queries:
        hs = "".join(
            f"<Hit><Hit_num>{i}</Hit_num><Hit_accession>{acc}</Hit_accession><Hit_def>{title}"
            f"</Hit_def><Hit_len>100</Hit_len><Hit_hsps><Hsp><Hsp_bit-score>70</Hsp_bit-score>"
            f"<Hsp_evalue>1e-20</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>40"
            f"</Hsp_query-to><Hsp_hit-from>50</Hsp_hit-from><Hsp_hit-to>11</Hsp_hit-to>"
            f"<Hsp_identity>{ident}</Hsp_identity><Hsp_align-len>40</Hsp_align-len></Hsp>"
            f"</Hit_hsps></Hit>" for i, (acc, title, ident) in enumerate(hits, 1))
        its.append(f"<Iteration><Iteration_query-def>{q}</Iteration_query-def>"
                   f"<Iteration_query-len>40</Iteration_query-len><Iteration_hits>{hs}"
                   f"</Iteration_hits></Iteration>")
    return ("<?xml version=\"1.0\"?><BlastOutput><BlastOutput_iterations>" + "".join(its)
            + "</BlastOutput_iterations></BlastOutput>")


def _genbank() -> str:
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq("ACGT" * 10), id="CP000001.1", name="CP000001",
                    description="Escherichia coli plasmid pX",
                    annotations={"molecule_type": "DNA"})
    rec.features.append(SeqFeature(FeatureLocation(0, 40), type="CDS",
                                   qualifiers={"gene": ["gyrA"],
                                               "product": ["DNA gyrase subunit A"]}))
    out = io.StringIO()
    from Bio import SeqIO
    SeqIO.write(rec, out, "genbank")
    return out.getvalue()


class FakeNcbi:
    """Answers like NCBI and keeps a clock that only moves when the client sleeps."""

    def __init__(self, xml):
        self.t, self.calls, self.xml, self.polls = 0.0, [], xml, 0

    def clock(self):
        return self.t

    def sleep(self, s):
        self.t += s

    def fetch(self, url, data):
        params = dict(urllib.parse.parse_qsl(data.decode() if data else url.split("?", 1)[1]))
        self.calls.append((self.t, params))
        if "efetch" in url:
            return _genbank()
        if params["CMD"] == "Put":
            return "<!--QBlastInfoBegin\n    RID = R1\n    RTOE = 15\nQBlastInfoEnd\n-->"
        if params.get("FORMAT_OBJECT") == "SearchInfo":
            self.polls += 1
            return "Status=READY\nThereAreHits=yes" if self.polls > 1 else "Status=WAITING"
        return self.xml if params["FORMAT_TYPE"] == "XML" else TEXT


def test_select_unitigs():
    rows = [("m1", 1, f"u{i:02d}", "A" * (20 + i)) for i in range(12)]
    rows += [("m2", 7, "u11", "A" * 31), ("m2", 8, "v1", "C" * 25)]
    m = pd.DataFrame(rows, columns=["model_id", "pattern_id", "unitig_id", "sequence"])
    sel = nc.select_unitigs(m, 10)
    assert sorted(sel["unitig_id"]) == sorted([f"u{i:02d}" for i in range(2, 12)] + ["v1"])


def test_client_keeps_ncbi_rules():
    fake = FakeNcbi(_xml([]))
    c = nc.Ncbi("me@example.org", fetch=fake.fetch, sleep=fake.sleep, clock=fake.clock)
    rid, xml, text = c.blast(">q\nACGT\n", {"WORD_SIZE": 11})
    assert rid == "R1" and text == TEXT
    times = [t for t, p in fake.calls]
    assert all(b - a >= 10 for a, b in zip(times, times[1:], strict=False))
    polls = [t for t, p in fake.calls if p.get("FORMAT_OBJECT") == "SearchInfo"]
    assert polls[1] - polls[0] >= 60
    put = fake.calls[0][1]
    assert put["PROGRAM"] == "blastn" and put["DATABASE"] == "nt" and put["EMAIL"] == "me@example.org"
    with pytest.raises(ValueError):
        nc.Ncbi("")


def test_parsing():
    hits = nc.parse_xml(_xml([("KU1", [("CP000001", "E. coli plasmid pX", 39),
                                        ("CP000002", "E. coli chromosome", 40)]), ("KU2", [])]))
    assert len(hits) == 2 and hits["identity"].tolist() == [97.5, 100.0]
    assert (hits["hit_start"].iloc[0], hits["hit_stop"].iloc[0]) == (11, 50)
    assert nc.nt_release(TEXT) == "Oct 1, 2026"
    assert nc.gene_and_product(_genbank()) == ("gyrA", "DNA gyrase subunit A")
    s = nc.summarize(hits, {("CP000001", 11, 50): ("gyrA", "DNA gyrase subunit A")})
    assert s.iloc[0]["plasmid_share"] == 0.5 and s.iloc[0]["gene"] == "gyrA"


def test_step_end_to_end(tmp_path, monkeypatch):
    seq = "ACGTTGCAAGGCTTACCGGATCCATGGTACCAGTACGTTGCA"
    uid = kb.unitig_id(seq)[0]
    config = {"paths_organism": {"panel_dir": str(tmp_path / "panel"),
                                 "card_layer_dir": str(tmp_path / "{organism}/{antibiotic}/card"),
                                 "context_dir": str(tmp_path / "{organism}/context")},
              "context": {"word_size": 11, "max_target_seqs": 50, "evalue": 10, "max_members": 10,
                          "batch_size": 100}}
    (tmp_path / "panel").mkdir()
    pd.DataFrame({"organism": ["ecoli"], "antibiotic": ["ciprofloxacin"],
                  "decision": ["included"]}).to_csv(tmp_path / "panel" / "panel_decisions.csv",
                                                    index=False)
    card = tmp_path / "ecoli" / "ciprofloxacin" / "card"
    card.mkdir(parents=True)
    pd.DataFrame({"pattern_id": [3], "sequence": [seq]}).to_csv(card / "card_unitigs.csv",
                                                               index=False)
    spec = importlib.util.spec_from_file_location("amrtest_18", PROJECT_ROOT / "scripts" /
                                                  "18_ncbi_context.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    fake = FakeNcbi(_xml([(uid, [("CP000001", "E. coli plasmid pX", 40)])]))
    client = nc.Ncbi("me@example.org", fetch=fake.fetch, sleep=fake.sleep, clock=fake.clock)
    out = tmp_path / "ecoli" / "context"
    r = m.query("ecoli", config, out, client)
    assert r == {"n_unitigs": 1, "n_searched_now": 1, "n_genbank_now": 1}
    put = fake.calls[0][1]
    assert put["ENTREZ_QUERY"] == "txid562[Organism:exp]" and put["HITLIST_SIZE"] == "50"
    assert m.entrez_query([550, 61645]) == "(txid550[Organism:exp] OR txid61645[Organism:exp])"
    n_calls = len(fake.calls)
    assert m.query("ecoli", config, out, client)["n_searched_now"] == 0   # from the cache
    assert len(fake.calls) == n_calls
    t = m.build("ecoli", out)
    from lib import contract  # the table follows the output contract
    assert contract.validate_csv(out / "unitig_context.csv",
                                 contract.load()["tables"]["unitig_context"]) == []
    assert json.loads((out / "context_summary.json").read_text())["n_unitigs"] == len(t)
    row = t.iloc[0]
    assert row["unitig_id"] == uid and row["nt_release"] == "Oct 1, 2026"
    assert row["gene"] == "gyrA" and row["plasmid_share"] == 1.0 and row["source"] == nc.SOURCE

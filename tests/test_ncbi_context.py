#!/usr/bin/env python3
"""NCBI context (lib.ncbi_context, 18_ncbi_context.py) with a fake NCBI: no network."""
import http.client
import importlib.util
import io
import json
import sys
import urllib.error
import urllib.parse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import knowledge_base as kb  # noqa: E402
from lib import ncbi_context as nc  # noqa: E402

pytestmark = pytest.mark.unit

SATURDAY = datetime(2026, 11, 7, 17, 0, tzinfo=timezone.utc)   # 12:00 US Eastern
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

    def __init__(self, xml, start=SATURDAY):
        self.t, self.calls, self.xml, self.polls, self.start = 0.0, [], xml, 0, start
        self.queries = []

    def clock(self):
        return self.t

    def now(self):
        return datetime.fromtimestamp(self.start.timestamp() + self.t, timezone.utc)

    def sleep(self, s):
        self.t += s

    def fetch(self, url, data):
        params = dict(urllib.parse.parse_qsl(data.decode() if data else url.split("?", 1)[1]))
        self.calls.append((self.t, params))
        if "efetch" in url:
            return _genbank()
        if params["CMD"] == "Put":
            self.queries.append(params["QUERY"])
            return "<!--QBlastInfoBegin\n    RID = R1\n    RTOE = 15\nQBlastInfoEnd\n-->"
        if params.get("FORMAT_OBJECT") == "SearchInfo":
            self.polls += 1
            return "Status=READY\nThereAreHits=yes" if self.polls > 1 else "Status=WAITING"
        return self.xml if params["FORMAT_TYPE"] == "XML" else TEXT


def client(fake, **kw):
    return nc.Ncbi("me@example.org", fetch=fake.fetch, sleep=fake.sleep, clock=fake.clock,
                   now=fake.now, **kw)


def test_select_unitigs():
    rows = [("m1", 1, f"u{i:02d}", "A" * (20 + i)) for i in range(12)]
    rows += [("m2", 7, "u11", "A" * 31), ("m2", 8, "v1", "C" * 25)]
    m = pd.DataFrame(rows, columns=["model_id", "pattern_id", "unitig_id", "sequence"])
    sel = nc.select_unitigs(m, 10)
    assert sorted(sel["unitig_id"]) == sorted([f"u{i:02d}" for i in range(2, 12)] + ["v1"])


def test_client_keeps_ncbi_rules():
    fake = FakeNcbi(_xml([]))
    c = client(fake)
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


def test_new_searches_wait_for_ncbi_hours():
    et = nc.EASTERN
    assert nc.until_off_hours(datetime(2026, 11, 9, 10, 0, tzinfo=et)) == 11 * 3600    # Monday
    assert nc.until_off_hours(datetime(2026, 11, 13, 20, 30, tzinfo=et)) == 1800       # Friday
    for t in [datetime(2026, 11, 9, 21, 0, tzinfo=et), datetime(2026, 11, 10, 4, 59, tzinfo=et),
              SATURDAY, datetime(2026, 11, 8, 23, 0, tzinfo=et)]:
        assert nc.until_off_hours(t) == 0
    assert nc.until_off_hours(datetime(2026, 11, 9, 15, 0, tzinfo=timezone.utc)) == 11 * 3600
    fake = FakeNcbi(_xml([]), start=datetime(2026, 11, 9, 10, 0, tzinfo=et))
    client(fake).blast(">q\nACGT\n", {})
    assert fake.calls[0][1]["CMD"] == "Put" and fake.calls[0][0] >= 11 * 3600
    fake = FakeNcbi(_xml([]), start=datetime(2026, 11, 9, 10, 0, tzinfo=et))
    client(fake, off_hours_only=False).blast(">q\nACGT\n", {})
    assert fake.calls[0][0] == 0
    # a search submitted just before 05:00 is followed to its end, not left for the night
    fake = FakeNcbi(_xml([]), start=datetime(2026, 11, 10, 4, 59, 50, tzinfo=et))
    client(fake).blast(">q\nACGT\n", {})
    assert fake.calls[-1][0] < 600


class Flaky(FakeNcbi):
    """Drops the connection the first ``n`` times the XML report is asked for."""

    def __init__(self, xml, n, error=None):
        super().__init__(xml)
        self.n, self.error = n, error

    def fetch(self, url, data):
        params = dict(urllib.parse.parse_qsl(data.decode() if data else url.split("?", 1)[1]))
        if params.get("FORMAT_TYPE") == "XML" and self.n > 0:
            self.n -= 1
            self.calls.append((self.t, params))
            raise self.error or http.client.RemoteDisconnected("closed")
        return super().fetch(url, data)


def test_client_sends_again_after_a_dropped_connection(capsys):
    fake = Flaky(_xml([]), 2)
    c = client(fake)
    rid, xml, text = c.blast(">q\nACGT\n", {"WORD_SIZE": 11})
    assert rid == "R1" and xml == _xml([]) and text == TEXT
    gets = [t for t, p in fake.calls if p.get("FORMAT_TYPE") == "XML"]
    assert len(gets) == 3 and gets[1] - gets[0] >= 60 and gets[2] - gets[1] >= 120
    assert sum(p.get("CMD") == "Put" for t, p in fake.calls) == 1      # the search is kept
    log = capsys.readouterr().err
    assert "NCBI blast XML report R1: RemoteDisconnected" in log and "attempt 2 of 11" in log
    assert "BLAST R1 submitted" in log
    fake = Flaky(_xml([]), 7)
    c = client(fake)
    c.blast(">q\nACGT\n", {})
    gets = [t for t, p in fake.calls if p.get("FORMAT_TYPE") == "XML"]
    assert len(gets) == 8 and max(b - a for a, b in zip(gets, gets[1:], strict=False)) <= 900


def test_client_gives_up():
    c = client(Flaky(_xml([]), 9), retries=2)
    with pytest.raises(http.client.RemoteDisconnected):
        c.blast(">q\nACGT\n", {})
    bad = urllib.error.HTTPError("u", 400, "Bad Request", None, None)  # type: ignore[arg-type]
    fake = Flaky(_xml([]), 1, bad)
    c = client(fake)
    with pytest.raises(urllib.error.HTTPError):
        c.blast(">q\nACGT\n", {})
    assert sum(p.get("FORMAT_TYPE") == "XML" for t, p in fake.calls) == 1
    assert nc.transient(urllib.error.HTTPError("u", 503, "x", None, None))  # type: ignore[arg-type]
    fake = FakeNcbi(_xml([]))
    fake.fetch = lambda url, data: "<html><p>Too many searches</p></html>"  # type: ignore[method-assign]
    with pytest.raises(nc.SearchFailed, match="Too many searches"):
        client(fake).blast(">q\nACGT\n", {})



def test_a_search_left_waiting_fails_after_the_limit():
    fake = FakeNcbi(_xml([]))
    fake.polls = -10 ** 6                       # NCBI answers WAITING for good
    with pytest.raises(nc.SearchFailed, match="still WAITING after 2.0 h"):
        client(fake, max_search_wait=2 * 3600).blast(">q\nACGT\n", {})
    assert 2 * 3600 <= fake.t < 2 * 3600 + 120
    assert not any(p.get("FORMAT_TYPE") == "XML" for t, p in fake.calls)


def test_parsing():
    hits = nc.parse_xml(_xml([("KU1", [("CP000001", "E. coli plasmid pX", 39),
                                        ("CP000002", "E. coli chromosome", 40)]), ("KU2", [])]))
    assert len(hits) == 2 and hits["identity"].tolist() == [97.5, 100.0]
    assert (hits["hit_start"].iloc[0], hits["hit_stop"].iloc[0]) == (11, 50)
    assert nc.nt_release(TEXT) == "Oct 1, 2026"
    assert nc.gene_and_product(_genbank()) == ("gyrA", "DNA gyrase subunit A")
    s = nc.summarize(hits, {("CP000001", 11, 50): ("gyrA", "DNA gyrase subunit A")})
    assert s.iloc[0]["plasmid_share"] == 0.5 and s.iloc[0]["gene"] == "gyrA"


SEQ = "ACGTTGCAAGGCTTACCGGATCCATGGTACCAGTACGTTGCA"


def _step(tmp_path, seqs, batch_size=100):
    """Step 18, and a config with one E. coli model whose CARD layer has ``seqs``."""
    config = {"paths_organism": {"panel_dir": str(tmp_path / "panel"),
                                 "card_layer_dir": str(tmp_path / "{organism}/{antibiotic}/card"),
                                 "context_dir": str(tmp_path / "{organism}/context")},
              "context": {"word_size": 11, "max_target_seqs": 50, "evalue": 10, "max_members": 10,
                          "batch_size": batch_size}}
    (tmp_path / "panel").mkdir(parents=True)
    pd.DataFrame({"organism": ["ecoli"], "antibiotic": ["ciprofloxacin"],
                  "decision": ["included"]}).to_csv(tmp_path / "panel" / "panel_decisions.csv",
                                                    index=False)
    card = tmp_path / "ecoli" / "ciprofloxacin" / "card"
    card.mkdir(parents=True)
    pd.DataFrame({"pattern_id": range(3, 3 + len(seqs)), "sequence": seqs}).to_csv(
        card / "card_unitigs.csv", index=False)
    spec = importlib.util.spec_from_file_location("amrtest_18", PROJECT_ROOT / "scripts" /
                                                  "18_ncbi_context.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m, config, tmp_path / "ecoli" / "context"


class FailingSearch(FakeNcbi):
    """Drops the connection every time a search for ``unitig`` is submitted, ``n`` times."""

    def __init__(self, xml, unitig, n):
        super().__init__(xml)
        self.unitig, self.n = unitig, n

    def fetch(self, url, data):
        params = dict(urllib.parse.parse_qsl(data.decode() if data else url.split("?", 1)[1]))
        if params.get("CMD") == "Put" and self.unitig in params["QUERY"] and self.n > 0:
            self.n -= 1
            self.queries.append(params["QUERY"])
            raise http.client.RemoteDisconnected("closed")
        return super().fetch(url, data)


def test_step_leaves_a_failing_batch_for_the_end(tmp_path):
    seqs = [SEQ, SEQ[::-1]]
    m, config, out = _step(tmp_path, seqs, batch_size=1)
    a, b = sorted(kb.unitig_id(s)[0] for s in seqs)
    fake = FailingSearch(_xml([(a, [("CP000001", "E. coli plasmid pX", 40)]),
                               (b, [("CP000001", "E. coli plasmid pX", 40)])]), a, 2)
    r = m.query("ecoli", config, out, client(fake, retries=1))
    assert r["n_searched_now"] == 2
    assert [b in q for q in fake.queries] == [False, False, True, False]   # a, a, b, then a
    m, config, out = _step(tmp_path / "again", seqs, batch_size=1)
    fake = FailingSearch(fake.xml, a, 4)
    with pytest.raises(SystemExit, match="NCBI failed 1 unitig"):
        m.query("ecoli", config, out, client(fake, retries=1))
    assert len(list((out / "cache" / "blast").glob("*/meta.json"))) == 1    # b is kept
    assert len(list((out / "cache" / "entrez").glob("*.gb"))) == 1


def test_step_end_to_end(tmp_path, monkeypatch):
    seq = SEQ
    uid = kb.unitig_id(seq)[0]
    m, config, out = _step(tmp_path, [seq])
    fake = FakeNcbi(_xml([(uid, [("CP000001", "E. coli plasmid pX", 40)])]))
    c = client(fake)
    r = m.query("ecoli", config, out, c)
    assert r == {"n_unitigs": 1, "n_searched_now": 1, "n_genbank_now": 1}
    put = fake.calls[0][1]
    assert put["ENTREZ_QUERY"] == "txid562[Organism:exp]" and put["HITLIST_SIZE"] == "50"
    assert m.entrez_query([550, 61645]) == "(txid550[Organism:exp] OR txid61645[Organism:exp])"
    n_calls = len(fake.calls)
    assert m.query("ecoli", config, out, c)["n_searched_now"] == 0   # from the cache
    assert len(fake.calls) == n_calls
    t = m.build("ecoli", out)
    from lib import contract  # the table follows the output contract
    assert contract.validate_csv(out / "unitig_context.csv",
                                 contract.load()["tables"]["unitig_context"]) == []
    assert json.loads((out / "context_summary.json").read_text())["n_unitigs"] == len(t)
    row = t.iloc[0]
    assert row["unitig_id"] == uid and row["nt_release"] == "Oct 1, 2026"
    assert row["gene"] == "gyrA" and row["plasmid_share"] == 1.0 and row["source"] == nc.SOURCE


class _Done:
    def __init__(self, rc, out=b"", err=b""):
        self.returncode, self.stdout, self.stderr = rc, out, err


def test_curl_returns_the_body_and_sends_post_data():
    seen = {}

    def run(cmd, input, capture_output, timeout):
        seen["cmd"], seen["input"] = cmd, input
        return _Done(0, b"RID = ABC" + nc.CURL_MARK.encode() + b"200")
    assert nc.curl("https://x/Blast.cgi", b"CMD=Put", run=run) == "RID = ABC"
    assert seen["cmd"][:3] == ["curl", "--data-binary", "@-"] and seen["input"] == b"CMD=Put"
    assert "-A" in seen["cmd"] and seen["cmd"][seen["cmd"].index("-A") + 1] == nc.TOOL


def test_curl_failures_are_judged_like_urllib():
    with pytest.raises(ConnectionError) as e:                 # 52: empty reply, the 10.11 failure
        nc.curl("https://x", None, run=lambda *a, **k: _Done(52, err=b"Empty reply from server"))
    assert nc.transient(e.value)
    with pytest.raises(urllib.error.HTTPError) as e:
        nc.curl("https://x", None, run=lambda *a, **k: _Done(0, b"busy" + nc.CURL_MARK.encode() + b"503"))
    assert e.value.code == 503 and nc.transient(e.value)
    with pytest.raises(urllib.error.HTTPError) as e:
        nc.curl("https://x", None, run=lambda *a, **k: _Done(0, b"no" + nc.CURL_MARK.encode() + b"400"))
    assert not nc.transient(e.value)

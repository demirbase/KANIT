"""NCBI Datasets v2: genome reports, taxonomy lineages and assembly FASTA files.

Transient failures (network, 429, 5xx) are retried with a growing wait, and calls are
spaced so that the rate stays under NCBI's limit without an API key (5 per second).
Used by the external validation (step 19) to download the external isolates' GenBank
assemblies and to check their NCBI taxon.
"""
from __future__ import annotations

import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from collections.abc import Callable, Iterable

from lib.bvbrc import USER_AGENT, _ssl_context

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
Fetch = Callable[[str], bytes]
BATCH = 50                                      # accessions per report request


def http_get(url: str, *, timeout: int = 300) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout, context=_ssl_context()) as r:
        return r.read()


class Datasets:
    def __init__(self, fetch: Fetch = http_get, *, tries: int = 5, wait: float = 5.0,
                 pause: float = 0.25, sleep: Callable[[float], None] = time.sleep):
        self.fetch, self.tries, self.wait, self.pause, self.sleep = fetch, tries, wait, pause, sleep

    def _get(self, path: str) -> bytes:
        for attempt in range(self.tries):
            try:
                body = self.fetch(f"{API}/{path}")
                self.sleep(self.pause)
                return body
            except urllib.error.HTTPError as e:
                if (e.code < 500 and e.code != 429) or attempt == self.tries - 1:
                    raise
            except (urllib.error.URLError, TimeoutError, ConnectionError):
                if attempt == self.tries - 1:
                    raise
            self.sleep(self.wait * 2 ** attempt)
        raise AssertionError("unreachable")

    def genome_reports(self, accessions: Iterable[str]) -> dict[str, dict]:
        """{accession: dataset report} of the assemblies NCBI knows; an accession it does
        not know (suppressed, withdrawn, mistyped) is absent."""
        accs = sorted({a.strip() for a in accessions if a and a.strip()})
        out: dict[str, dict] = {}
        for i in range(0, len(accs), BATCH):
            batch = ",".join(urllib.parse.quote(a) for a in accs[i:i + BATCH])
            doc = json.loads(self._get(f"genome/accession/{batch}/dataset_report?page_size=1000"))
            for r in doc.get("reports", []):
                out[r["accession"]] = r
        return out

    def lineages(self, tax_ids: Iterable[int]) -> dict[int, set[int]]:
        """{tax id: the taxon and all its ancestors}."""
        ids = sorted({int(t) for t in tax_ids})
        out: dict[int, set[int]] = {}
        for i in range(0, len(ids), BATCH):
            batch = ",".join(str(t) for t in ids[i:i + BATCH])
            doc = json.loads(self._get(f"taxonomy/taxon/{batch}/dataset_report"))
            for r in doc.get("reports", []):
                t = r["taxonomy"]
                out[int(t["tax_id"])] = {int(t["tax_id"]), *(int(p) for p in t.get("parents", []))}
        return out

    def fasta(self, accession: str) -> bytes:
        """The genomic FASTA of one assembly (the GENOME_FASTA of its package)."""
        body = self._get(f"genome/accession/{urllib.parse.quote(accession)}/download"
                         "?include_annotation_type=GENOME_FASTA")
        with zipfile.ZipFile(io.BytesIO(body)) as z:
            names = [n for n in z.namelist() if n.endswith("_genomic.fna")]
            if len(names) != 1:
                raise ValueError(f"{accession}: {len(names)} genomic FASTA files in the package")
            return z.read(names[0])


def report_numbers(report: dict) -> tuple[str, str, int | None, str]:
    """(contigs, total length, tax id, organism name) of a genome report, as strings for
    bvbrc.check_fasta ('' when absent)."""
    stats = report.get("assembly_stats", {})
    org = report.get("organism", {})
    tax = org.get("tax_id")
    return (str(stats.get("number_of_contigs", "")), str(stats.get("total_sequence_length", "")),
            int(tax) if tax is not None else None, str(org.get("organism_name", "")))

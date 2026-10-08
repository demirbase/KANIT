#!/usr/bin/env python3
"""Outer folds and permutation chunks per SLURM job (scripts/lib/packing.py)."""
import math
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from lib import packing  # noqa: E402
from lib.config import load_config  # noqa: E402

pytestmark = pytest.mark.unit


def test_units_per_job():
    cfg = load_config()["packing"]
    # the E. coli probe (5,681 genomes, 4.70 M patterns): two units fill a job
    assert packing.unit_memory_gb(5681, 4_700_000, cfg) > 60
    assert packing.per_job(5681, 4_700_000, cfg) == 2
    # small models: as many as the memory holds (a fixed 4 GB each keeps them under a
    # hamsi node's 56)
    assert packing.per_job(1384, 500_000, cfg) == 29
    assert packing.per_job(1204, 1_118_452, cfg) == 22                   # pilot, imipenem
    assert packing.per_job(100, 1_000, {**cfg, "max_per_job": 20}) == 20   # the cap
    # the estimate covers the peaks measured in the pilot (A. baumannii, 2026-10-06; genomes,
    # patterns, GB per unit)
    for n, p, peak in [(968, 894_710, 5.8), (760, 759_393, 5.4), (435, 614_743, 3.7),
                       (784, 746_163, 4.6), (689, 758_695, 4.3), (750, 790_436, 4.1)]:
        assert packing.unit_memory_gb(n, p, cfg) >= peak
    # in between: as many as the memory holds
    mid = {**cfg, "max_per_job": 99}
    k = packing.per_job(3000, 2_000_000, mid)
    assert k * packing.unit_memory_gb(3000, 2_000_000, mid) <= mid["budget_gb"] < \
        (k + 1) * packing.unit_memory_gb(3000, 2_000_000, mid)
    assert packing.per_job(10**6, 10**8, cfg) == 1                     # never below one


def test_job_memory_holds_the_budget():
    """The units of a job fit in the memory of the labels that run them on every partition
    that runs them (capped there, as in memoryFor, at 95% of a node), and on barbun that
    memory needs no more cores than a barbun job gets anyway (conf/truba.config)."""
    text = (ROOT / "conf" / "truba.config").read_text()
    parts = {name: [int(x) for x in rest] for name, *rest in re.findall(
        r"(\w+)\s*:\s*\[cpus: (\d+), min_cpus: (\d+), mb_per_cpu: (\d+)", text)}
    budget_gb = load_config()["packing"]["budget_gb"]
    for label in ("cv_unit", "permutation"):
        gib = int(re.search(rf"withLabel: {label} \{{\s*(?:queue .*\n\s*)?"
                            rf"memory = \{{ memoryFor\((\d+),", text).group(1))
        for name in ("barbun", "hamsi", "orfoz"):
            cpus, _, mb_per_cpu = parts[name]
            mb = min(gib * 1024, math.floor(cpus * mb_per_cpu * 0.95))
            assert budget_gb * 1e9 <= mb * 2**20, (label, name)
        cpus, min_cpus, mb_per_cpu = parts["barbun"]
        assert math.ceil(gib * 1024 / mb_per_cpu) <= min_cpus <= cpus, label

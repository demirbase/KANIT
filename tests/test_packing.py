#!/usr/bin/env python3
"""Outer folds and permutation chunks per SLURM job (scripts/lib/packing.py)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from lib import packing  # noqa: E402
from lib.config import load_config  # noqa: E402

pytestmark = pytest.mark.unit


def test_units_per_job():
    cfg = load_config()["packing"]
    # the E. coli probe (5,681 genomes, 4.70 M patterns): one unit fills a job
    assert packing.unit_memory_gb(5681, 4_700_000, cfg) > 60
    assert packing.per_job(5681, 4_700_000, cfg) == 1
    # a small model: as many as the cores allow
    assert packing.per_job(1384, 500_000, cfg) == cfg["max_per_job"]
    # in between: as many as the memory holds
    mid = {**cfg, "max_per_job": 99}
    k = packing.per_job(3000, 2_000_000, mid)
    assert k * packing.unit_memory_gb(3000, 2_000_000, mid) <= mid["budget_gb"] < \
        (k + 1) * packing.unit_memory_gb(3000, 2_000_000, mid)
    assert packing.per_job(10**6, 10**8, cfg) == 1                     # never below one

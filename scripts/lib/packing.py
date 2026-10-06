"""How many outer folds (CV units) or label permutation chunks of one model share a SLURM
job. Execution only: the units, their seeds and settings are unchanged, and every CV unit
records the number of threads it ran with.

TRUBA's barbun partition gives a job at least 20 cores, and a unit gains little from more than
one thread (pilot, A. baumannii, 2026-10-06: 242% CPU with 20 threads, 85-138% with two).
Several units therefore run side by side in one job, each with its share of the cores, as many
as the job's memory holds: a unit's peak memory grows with its training matrix
(genomes × patterns; E. coli probe 2026-10-02: 4.70 M patterns × 5,681 genomes, 76 GiB), and
the rest is about 3-4 GB (pilot, units in jobs of up to six: 3.1-5.8 GB each, 2.6-3.9 GB of it
besides the matrix).
"""
from __future__ import annotations

import math


def unit_memory_gb(n_genomes: int, n_patterns: int, cfg: dict) -> float:
    """Expected peak memory of one unit: a fixed part and the outer training set's matrix
    (four fifths of the model's genomes)."""
    cells = 0.8 * n_genomes * n_patterns
    return float(cfg["base_gb"]) + float(cfg["bytes_per_cell"]) * cells / 1e9


def per_job(n_genomes: int, n_patterns: int, cfg: dict) -> int:
    """Units of the model in one job: as many as ``budget_gb`` holds, at most ``max_per_job``,
    at least one."""
    k = math.floor(float(cfg["budget_gb"]) / unit_memory_gb(n_genomes, n_patterns, cfg))
    return max(1, min(int(cfg["max_per_job"]), k))

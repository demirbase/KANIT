"""How many outer folds (CV units) or label permutation chunks of one model share a SLURM
job. Execution only: the analysis of each unit is unchanged.

TRUBA's barbun partition gives a job at least 20 cores, and one unit keeps about 2.4 of them
busy (pilot, A. baumannii, 2026-10-06: 113 min, CPU 242%, peak 6.5 GB on average). Several
units therefore run side by side in one job, each with its share of the cores, as many as
the job's memory holds: a unit's peak memory grows with its training matrix (genomes ×
patterns; E. coli probe 2026-10-02: 4.70 M patterns × 5,681 genomes, 76 GiB).
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

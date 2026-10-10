"""How many outer folds (CV units) or label permutation chunks of one model share a SLURM
job. Execution only: the units, their seeds and settings are unchanged, and every CV unit
records the number of threads it ran with.

TRUBA's barbun partition gives a job at least 20 cores, and a unit keeps few of them busy
(pilot, A. baumannii, 2026-10-06: 242% CPU with 20 threads; 85-138% with two, in a median 1.2
to 1.6 times the time). Packed, a unit used about a quarter of the core hours, in jobs that
held one to six new units.
Several units therefore run side by side in one job, each with its share of the cores, as many
as the job's memory holds. A unit's peak memory grows with its training matrix (genomes ×
patterns; E. coli probe 2026-10-02: 4.70 M patterns × 5,681 genomes, 76 GiB) or, in a model of
fewer genomes, with its number of patterns (full run 2026-10-10, peaks per unit: P. aeruginosa
meropenem, 913 genomes × 1.62 M patterns, 8-20 GB; E. coli tobramycin, 961 × 2.18 M, 14-32 GB;
the pilot's A. baumannii units, under 1.1 M patterns, 3.1-5.8 GB). The estimate takes the larger
of the two, set near the units' mean peak: units rarely peak together, and a job that runs out
of memory is retried with fewer units at once (modules/local/common.nf, parallel).
"""
from __future__ import annotations

import math


def unit_memory_gb(n_genomes: int, n_patterns: int, cfg: dict) -> float:
    """Expected peak memory of one unit: a fixed part and the larger of the outer training
    set's matrix (four fifths of the model's genomes) and a part per pattern."""
    cells = 0.8 * n_genomes * n_patterns
    return float(cfg["base_gb"]) + max(float(cfg["bytes_per_cell"]) * cells / 1e9,
                                       float(cfg["gb_per_million_patterns"]) * n_patterns / 1e6)


def per_job(n_genomes: int, n_patterns: int, cfg: dict) -> int:
    """Units of the model in one job: as many as ``budget_gb`` holds, at most ``max_per_job``,
    at least one."""
    k = math.floor(float(cfg["budget_gb"]) / unit_memory_gb(n_genomes, n_patterns, cfg))
    return max(1, min(int(cfg["max_per_job"]), k))

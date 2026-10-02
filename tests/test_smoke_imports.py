#!/usr/bin/env python3
"""
Smoke tests — load every numbered script and assert it imports cleanly.

This is the cheapest, highest-value safety net: it catches syntax errors, bad
imports, broken config wiring, and unresolved paths in SECONDS — instead of
discovering them hours into a multi-day pipeline run. Scripts that need a heavy
dependency missing from this environment (xgboost/optuna) are skipped, so the
same suite runs fully in your real environment.
"""

import pytest

PIPELINE_SCRIPTS = [
    "00a_download_bvbrc.py",
    "00_prepare_metadata.py",
    "01_data_validation.py",
    "01b_data_validation.py",
    "02c_lineage_poppunk.py",
    "02d_genome_qc.py",
    "02e_panel.py",
    "03u_unitig_matrix.py",
    "04_optimization.py",
    "05_model_training.py",
    "06_evaluation.py",
    "07_explainability.py",
    "07b_feature_stability.py",
    "08_blast_annotation.py",
    "09_biological_summary.py",
    "10_unitig_background_frequency.py",
    "11_variant_snp_check.py",
]


@pytest.mark.smoke
@pytest.mark.parametrize("script", PIPELINE_SCRIPTS)
def test_script_imports(load_script, script):
    """Each script loads without syntax/import/config errors (or skips on missing dep)."""
    mod = load_script(script)
    assert mod is not None


@pytest.mark.smoke
def test_lib_package_imports():
    """The shared lib package and its public API import cleanly."""
    from lib import chunking, config, io_utils, registry, run_metadata  # noqa: F401
    assert callable(config.resolve_path)
    assert callable(registry.load_antibiotic_classes)
    assert callable(chunking.get_y_chunk)
    assert callable(io_utils.run_command)
    assert callable(run_metadata.make_run_id)


# ============================================================================
# Makefile — ML_AMR_Prediction_v2 developer & pipeline shortcuts
# ============================================================================
# `make help` lists targets. PYTHON can be overridden, e.g.:
#     make test PYTHON=/opt/anaconda3/envs/bitirme_vol2/bin/python
# ============================================================================
PYTHON ?= python
ORG    ?= ecoli
AB     ?= ampicillin

.DEFAULT_GOAL := help
.PHONY: help setup dev-install lint format typecheck test test-all \
        pipeline data train biology tables clean-pyc

# Artefact paths. Override on the command line if the KB lives elsewhere,
# e.g. `make tables KB=$$AMR_WORK/results/kb/amrk.db`.
KB      ?= results/kb/amrk.db
RESULTS ?= results
TABLES  ?= results/tables
FIGURES ?= results/figures

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
	  awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

setup:  ## Create the conda environment (external tools + Python deps)
	conda env create -f environment.yml

dev-install:  ## Install the dev/QA toolchain (ruff, mypy, pre-commit, pytest)
	$(PYTHON) -m pip install -e ".[dev]"
	pre-commit install

lint:  ## Lint with ruff
	ruff check .

format:  ## Auto-format with ruff
	ruff format .

typecheck:  ## Type-check the lib/ core with mypy
	mypy

test:  ## Run unit + smoke tests (fast)
	$(PYTHON) -m pytest -ra

test-all:  ## Run everything incl. the opt-in integration tests
	$(PYTHON) -m pytest -m "unit or smoke or integration" -ra

pipeline:  ## Run the analysis core (01->10) for ORG/AB from config
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --antibiotic $(AB)

data:  ## Acquire BV-BRC data + build the phenotype matrix (00a, 00)
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --only 00a 00

train:  ## HPO + train + evaluate (04, 05, 06)
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --only 04 05 06

biology:  ## Stability + explainability + BLAST + reports (07b..11)
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --only 07b 07 08 09 10 11

tables:  ## Rebuild the tidy tables from the KB (kb_tables + H3 + CV comparison)
	$(PYTHON) scripts/kb_tables.py --db $(KB) --results $(RESULTS) --out $(TABLES)
	$(PYTHON) scripts/17_h3_gene_family_overlap.py --db $(KB) --tables $(TABLES) --figures $(FIGURES)
	$(PYTHON) scripts/kb_cv_comparison.py --tables $(TABLES) --results $(RESULTS) --out $(FIGURES)
	$(PYTHON) scripts/18_novel_ncbi_context.py --kb $(KB) --results-root $(RESULTS) --out $(TABLES)
	$(PYTHON) scripts/kb_fair_mapping.py --db $(KB) --out $(TABLES)

clean-pyc:  ## Remove Python caches
	find . -type d -name __pycache__ -prune -exec rm -rf {} + ; \
	find . -type f -name '*.py[co]' -delete

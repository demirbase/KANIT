# KANIT

**K**nowledgebase for **A**ntimicrobial Resistance from machine-learned u**N**itigs with
**I**ntegrated evidence **T**iering.

KANIT is a knowledge base of genomic biomarkers associated with antimicrobial resistance in
the ESKAPEE pathogens. The biomarkers are unitigs, paths of a compacted de Bruijn graph built
from bacterial genomes, that gradient-boosted models select under lineage-aware
cross-validation. Each biomarker is examined by independent lines of evidence, which a fixed
rule combines into one evidence grade; a grade states association with resistance, not that
the biomarker causes it.

> **Status:** under development. The first release, v1.0.0, will include the knowledge base,
> the analysis workflow and the methods documentation.

## Repository layout

| Path | Contents |
|:--|:--|
| `scripts/` | the numbered analysis steps and the knowledge-base tools |
| `scripts/lib/` | shared code: configuration, registries, schema, I/O |
| `config/` | configuration, and the organism and antibiotic registries |
| `tests/` | unit and smoke tests |
| `*.def`, `environment*.yml` | container definitions and their environments |
| `slurm/` | HPC job scripts |

## Tests

```bash
pip install -e ".[dev]"     # or: conda env create -f environment.yml
pytest                       # unit and smoke tests
```

## Licence

The code is released under the MIT License (see `LICENSE`).

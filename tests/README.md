# Test Suite

A layered safety net so you can validate the pipeline in **seconds/minutes**
instead of rerunning the full multi-day job. None of the 9 numbered scripts are
modified — the suite is purely additive.

## Layers

| Layer | Marker | Speed | What it catches |
|-------|--------|-------|-----------------|
| **Smoke** | `smoke` | seconds | Syntax / import / config-wiring / unresolved-path breakage in every script |
| **Unit** | `unit` | seconds | Logic regressions in the bug-prone pure functions (R/S counting, BLAST tiers, √p colsample, bootstrap CI, threshold, chunking, registry, path resolution) |

## How to run

```bash
# Fast feedback — default. Runs unit + smoke, NEVER touches your config.
pytest

# Just one layer
pytest -m unit
pytest -m smoke
```

`pytest` defaults to `-m "not integration and not slow"` (see `pytest.ini`), so
heavy tests are **opt-in only**.

## Environments

In an environment without `xgboost`/`optuna` (e.g. CI lint box), the
xgboost-dependent unit tests **skip gracefully**; smoke + the remaining unit tests
still run. In your real training environment
everything runs.

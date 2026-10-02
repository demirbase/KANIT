#!/usr/bin/env python3
"""Step 04 — nested cross-validation and the final model of one panel pair.

Subcommands, so that a workflow can run every outer fold as its own task:

  folds    outer folds of both arms (lib.folds): folds.csv, seeds.csv,
           fold_composition.csv, cv_design.json
  unit     search, tree count, fit and prediction of one outer fold
           (--arm --repeat --fold): units/<arm>_r<repeat>_f<fold>/
  final    the final model on every genome: lineage-grouped inner splits, seed 0
  metrics  every unit of every evaluable split must be present; pooled
           out-of-fold metrics and lineage-cluster bootstrap intervals
  all      folds, every unit, final and metrics in one process

A model whose lineage-aware arm has a repeat without a valid seed is not
evaluable: `folds` records it and the other subcommands exit with status 3.
Outputs go to config paths_organism.cv_dir.

Usage:
    python scripts/04_nested_cv.py all --organism ecoli --antibiotic ampicillin
    python scripts/04_nested_cv.py unit --organism ecoli --antibiotic ampicillin \\
        --arm lineage_aware --repeat 1 --fold 0
"""
import argparse
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import folds, oof_metrics, train  # noqa: E402
from lib.config import get_target, load_config, resolve_path  # noqa: E402
from lib.matrix_store import ModelMatrix  # noqa: E402
from lib.run_metadata import peak_rss_gb  # noqa: E402

NOT_EVALUABLE = 3


def _now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def unit_name(arm, repeat, fold):
    return f"{arm}_r{repeat}_f{fold}"


def _data(mm):
    return (mm.genomes["Genome ID"].astype(str).to_numpy(), mm.labels.astype(int),
            mm.genomes["lineage"].astype(str).to_numpy())


# ---------------------------------------------------------------------------
def run_folds(mm, out_dir, cv):
    ids, y, groups = _data(mm)
    splits = folds.assign_all(y, groups, n_repeats=cv["n_repeats"], n_folds=cv["n_folds"],
                              min_minority=cv["min_minority_per_test_fold"],
                              max_attempts=cv["max_seed_attempts"])
    out_dir.mkdir(parents=True, exist_ok=True)
    folds.folds_table(ids, splits).to_csv(out_dir / "folds.csv", index=False)
    seeds = folds.seeds_table(splits)
    seeds.to_csv(out_dir / "seeds.csv", index=False)
    folds.composition_table(y, groups, splits, cv["n_folds"]).to_csv(
        out_dir / "fold_composition.csv", index=False)
    import sklearn
    share = folds.largest_lineage_share(groups)
    aware = seeds[seeds["arm"] == folds.LINEAGE_AWARE]
    design = {
        "created_at": _now(), "n_genomes": int(len(y)), "n_resistant": int(y.sum()),
        "n_susceptible": int((y == 0).sum()), "minority_label": folds.minority_label(y),
        "n_lineages": int(len(set(groups.tolist()))), "largest_lineage_share": share,
        "flag_largest_lineage": bool(share > 1 / cv["n_folds"]),
        "evaluable": bool(aware["evaluable"].all()),
        "lineage_blind_evaluable": bool(seeds.loc[seeds["arm"] == folds.LINEAGE_BLIND,
                                                  "evaluable"].all()),
        "cv": cv, "scikit_learn_version": sklearn.__version__,
    }
    _write_json(out_dir / "cv_design.json", design)
    return design


def _design(out_dir):
    path = out_dir / "cv_design.json"
    if not path.exists():
        sys.exit(f"ERROR: {path} not found; run `folds` first.")
    design = json.loads(path.read_text())
    if not design["evaluable"]:
        print("Model not evaluable: a lineage-aware repeat has no seed meeting the "
              "class-balance rule (see seeds.csv).")
        sys.exit(NOT_EVALUABLE)
    return design


def _split(out_dir, arm, repeat, ids):
    ft = pd.read_csv(out_dir / "folds.csv", dtype={"genome_id": str})
    ft = ft[(ft["arm"] == arm) & (ft["repeat"] == repeat)]
    if ft.empty:
        sys.exit(f"ERROR: no folds for {arm} repeat {repeat} (not evaluable?).")
    fold_of = pd.Series(ft["fold"].to_numpy(), index=ft["genome_id"]).reindex(ids)
    if fold_of.isna().any():
        sys.exit("ERROR: folds.csv does not cover the model's genomes.")
    seeds = pd.read_csv(out_dir / "seeds.csv")
    seed = seeds.loc[(seeds["arm"] == arm) & (seeds["repeat"] == repeat), "seed"].iloc[0]
    return fold_of.to_numpy(dtype=int), int(seed)


def _save(result, unit_dir, extra):
    tmp = unit_dir.with_name(unit_dir.name + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    result["booster"].save_model(tmp / "model.ubj")
    result["trials"].to_csv(tmp / "trials.csv", index=False)
    for name, df in extra.pop("tables", {}).items():
        df.to_csv(tmp / name, index=False)
    _write_json(tmp / "record.json", {**result["record"], **extra})
    shutil.rmtree(unit_dir, ignore_errors=True)
    tmp.rename(unit_dir)


def run_unit(mm, out_dir, arm, repeat, fold, hpo, threads):
    _design(out_dir)
    ids, y, groups = _data(mm)
    fold_of, seed = _split(out_dir, arm, repeat, ids)
    tr, te = np.flatnonzero(fold_of != fold), np.flatnonzero(fold_of == fold)
    t0 = time.time()
    result = train.train(mm, tr, y, groups, arm, seed, hpo, threads)
    p = train.predict(result["booster"], mm, te)
    oof = pd.DataFrame({"genome_id": ids[te], "y": y[te], "p": p})
    _save(result, out_dir / "units" / unit_name(arm, repeat, fold),
          {"repeat": repeat, "fold": fold, "n_test": int(len(te)), "finished_at": _now(),
           "seconds": round(time.time() - t0, 1), "peak_rss_gb": round(peak_rss_gb(), 2),
           "tables": {"oof.csv": oof}})


def run_final(mm, out_dir, hpo, threads):
    _design(out_dir)
    _, y, groups = _data(mm)
    t0 = time.time()
    result = train.train(mm, np.arange(len(y)), y, groups, folds.LINEAGE_AWARE, 0, hpo, threads)
    _save(result, out_dir / "final",
          {"finished_at": _now(), "seconds": round(time.time() - t0, 1),
           "peak_rss_gb": round(peak_rss_gb(), 2)})


def run_metrics(mm, out_dir, cv):
    _design(out_dir)
    ids, _, groups = _data(mm)
    seeds = pd.read_csv(out_dir / "seeds.csv")
    expected = [(s.arm, int(s.repeat), k) for s in seeds.itertuples() if s.evaluable
                for k in range(cv["n_folds"])]
    missing = [unit_name(*u) for u in expected
               if not (out_dir / "units" / unit_name(*u) / "record.json").exists()]
    if missing:
        sys.exit(f"ERROR: {len(missing)} outer fold(s) not finished: {', '.join(missing[:10])}")
    parts = []
    for arm, repeat, fold in expected:
        d = pd.read_csv(out_dir / "units" / unit_name(arm, repeat, fold) / "oof.csv",
                        dtype={"genome_id": str})
        parts.append(d.assign(arm=arm, repeat=repeat, fold=fold))
    oof = pd.concat(parts, ignore_index=True)[["arm", "repeat", "fold", "genome_id", "y", "p"]]
    counts = oof.groupby(["arm", "repeat"])["genome_id"].agg(["size", "nunique"])
    if not ((counts["size"] == len(ids)) & (counts["nunique"] == len(ids))).all():
        sys.exit("ERROR: an arm/repeat does not predict every genome exactly once.")
    oof.to_csv(out_dir / "oof_predictions.csv", index=False)

    by_arm = {arm: d.drop(columns="arm") for arm, d in oof.groupby("arm")}
    lineage_of = dict(zip(ids, groups, strict=True))
    diff = ((folds.LINEAGE_BLIND, folds.LINEAGE_AWARE)
            if len(by_arm) == 2 else None)
    ci = oof_metrics.cluster_bootstrap(by_arm, lineage_of, n_boot=cv["n_bootstrap"],
                                       seed=0, diff=diff)
    tables = {"repeat_metrics.csv": oof_metrics.per_repeat,
              "fold_metrics.csv": oof_metrics.per_fold}
    for name, fn in tables.items():
        pd.concat([fn(d, cv["threshold"]) if fn is oof_metrics.per_repeat else fn(d)
                   for d in by_arm.values()], keys=list(by_arm), names=["arm", None]
                  ).reset_index(level=0).to_csv(out_dir / name, index=False)
    pd.concat([oof_metrics.reliability(d, cv["reliability_bins"]) for d in by_arm.values()],
              keys=list(by_arm), names=["arm", None]).reset_index(level=0).to_csv(
        out_dir / "reliability.csv", index=False)
    metrics = {"created_at": _now(),
               "arms": {arm: {**oof_metrics.summary(d, cv["threshold"]), "ci": ci[arm]}
                        for arm, d in by_arm.items()},
               "bootstrap": {k: v for k, v in ci.items() if k not in by_arm}}
    _write_json(out_dir / "metrics.json", metrics)
    return metrics


# ---------------------------------------------------------------------------
def main():
    config = load_config()
    default_org, default_ab = get_target(config=config)
    ap = argparse.ArgumentParser(description="Nested cross-validation of one panel pair.")
    ap.add_argument("command", choices=["folds", "unit", "final", "metrics", "all"])
    ap.add_argument("--organism", default=default_org)
    ap.add_argument("--antibiotic", default=default_ab)
    ap.add_argument("--arm", choices=folds.ARMS)
    ap.add_argument("--repeat", type=int)
    ap.add_argument("--fold", type=int)
    ap.add_argument("--threads", type=int, default=int(config["preprocessing"]["threads"]))
    args = ap.parse_args()

    cv, hpo = config["cv"], config["hpo"]
    mm = ModelMatrix(resolve_path("matrix_dir", organism=args.organism,
                                  antibiotic=args.antibiotic, config=config))
    out_dir = resolve_path("cv_dir", organism=args.organism, antibiotic=args.antibiotic,
                           config=config)
    print(f"NESTED CV — {args.organism} / {args.antibiotic} — {args.command}")
    if args.command in ("folds", "all"):
        d = run_folds(mm, out_dir, cv)
        print(f"  evaluable: {d['evaluable']}; largest lineage {d['largest_lineage_share']:.1%}")
    if args.command == "unit":
        if args.arm is None or args.repeat is None or args.fold is None:
            sys.exit("ERROR: `unit` needs --arm, --repeat and --fold.")
        run_unit(mm, out_dir, args.arm, args.repeat, args.fold, hpo, args.threads)
    if args.command == "all":
        seeds = pd.read_csv(out_dir / "seeds.csv")
        for s in seeds.itertuples():
            if s.evaluable:
                for k in range(cv["n_folds"]):
                    run_unit(mm, out_dir, s.arm, int(s.repeat), k, hpo, args.threads)
                    print(f"  ✓ {unit_name(s.arm, int(s.repeat), k)}")
    if args.command in ("final", "all"):
        run_final(mm, out_dir, hpo, args.threads)
        print("  ✓ final model")
    if args.command in ("metrics", "all"):
        m = run_metrics(mm, out_dir, cv)
        for arm, s in m["arms"].items():
            print(f"  {arm}: ROC-AUC {s['roc_auc']:.3f} "
                  f"[{s['ci']['low']:.3f}, {s['ci']['high']:.3f}]")


if __name__ == "__main__":
    main()

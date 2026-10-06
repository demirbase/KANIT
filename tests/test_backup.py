#!/usr/bin/env python3
"""Verified backup (scripts/backup.py, scripts/backup_upload.sh) against a stand-in rclone."""
import csv
import importlib.util
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
pytestmark = pytest.mark.unit

FAKE_RCLONE = r"""#!/usr/bin/env bash
cmd=$1; shift
case $cmd in
  copy)  src=$1; dst=$2; mkdir -p "$dst"; n=$(basename "$src")
         if [ "${FAKE_FAIL_ONCE:-}" = "$n" ] && [ ! -e "$dst/.failed_$n" ]; then
             touch "$dst/.failed_$n"; exit 1; fi
         if [ "${FAKE_HANG_ONCE:-}" = "$n" ] && [ ! -e "$dst/.hung_$n" ]; then
             touch "$dst/.hung_$n"; exec sleep 30; fi
         if [ "${FAKE_CORRUPT:-}" = "$n" ]; then echo corrupted > "$dst/$n"; else cp "$src" "$dst/"; fi ;;
  check) src=$1; dst=$2; shift 2; inc=""
         while [ $# -gt 0 ]; do [ "$1" = --include-from ] && inc=$2; shift; done
         while read -r f; do cmp -s "$src/$f" "$dst/$f" || { echo "differ: $f" >&2; exit 1; }; done < "$inc" ;;
  lsf)   ls -1 "$1"; [ -n "${FAKE_DUP:-}" ] && echo "$FAKE_DUP"; true ;;
esac
"""


def _script():
    spec = importlib.util.spec_from_file_location("amrtest_backup",
                                                  PROJECT_ROOT / "scripts" / "backup.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m                     # dataclasses look their module up there
    spec.loader.exec_module(m)
    return m


def _tree(root: Path):
    files = {"results/ecoli/ampicillin/cv/oof.csv": "a\n1\n",
             "results/ecoli/ampicillin/grades/g.csv": "g\n",
             "results/panel/panel_decisions.csv": "p\n",
             "results/readme.txt": "r\n",
             "data/processed/ecoli/unitig_store/meta.json": "{}",
             "data/processed/ecoli/unitig_store/call/unitigs.rtab": "huge\n",
             "runs/nextflow/20261005T100000/run_manifest.json": "{}",
             "runs/nextflow/20261006T100000/run_manifest.json": "{}"}
    for rel, text in files.items():
        f = root / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(text)


def test_units_and_pack(tmp_path):
    m = _script()
    _tree(tmp_path)
    units = {(u.path, u.mode): u for u in m.units(tmp_path, current_run="20261006T100000")}
    assert set(units) == {("results", "files"), ("results/ecoli/ampicillin", "tree"),
                          ("results/panel", "files"),
                          ("data/processed/ecoli/unitig_store", "tree"),
                          ("runs/nextflow/20261005T100000", "tree")}
    assert units[("data/processed/ecoli/unitig_store", "tree")].names() == [
        "data/processed/ecoli/unitig_store/meta.json"]        # the unitig-caller call is left out
    stage, ledger = tmp_path / "stage", tmp_path / "ledger.tsv"
    rows = m.pack(stage, ledger, current_run="20261006T100000", threads=1, project=tmp_path)
    assert len(rows) == 5
    a = next(r for r in rows if r["unit"] == "results/ecoli/ampicillin")
    with tarfile.open(stage / a["archive"], "r:gz") as t:
        assert sorted(x.name for x in t.getmembers() if x.isfile()) == [
            "results/ecoli/ampicillin/cv/oof.csv", "results/ecoli/ampicillin/grades/g.csv"]
    assert a["archive"] == "results__ecoli__ampicillin.tar.gz" and a["n_files"] == 2
    # a unit already verified with the same fingerprint is not packed again
    with open(ledger, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=m.LEDGER_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows({**r, "remote": "x", "run": "r1", "verified_at": "2026-10-05T10:00:00Z"}
                    for r in rows)
    (tmp_path / "results/ecoli/ampicillin/cv/oof.csv").write_text("a\n2\n")
    again = m.pack(stage, ledger, current_run="20261006T100000", threads=1, project=tmp_path)
    assert [r["unit"] for r in again] == ["results/ecoli/ampicillin"]


def test_links_are_archived_as_links(tmp_path):
    """The genome QC links every assembly into one folder: the archive holds the links,
    and the assemblies are backed up where they lie."""
    m = _script()
    genome = tmp_path / "data/raw/ecoli/genomes/562.1.fna"
    genome.parent.mkdir(parents=True)
    genome.write_text(">c\nACGT\n")
    qc = tmp_path / "results/ecoli/global_exploration/genome_qc"
    (qc / "inputs").mkdir(parents=True)
    (qc / "inputs/562.1.fna").symlink_to(genome.resolve())
    (qc / "report.csv").write_text("r\n")
    rows = m.pack(tmp_path / "stage", tmp_path / "ledger.tsv", current_run=None, threads=1,
                  project=tmp_path)
    row = next(r for r in rows if r["unit"] == "results/ecoli/global_exploration")
    with tarfile.open(tmp_path / "stage" / row["archive"], "r:gz") as t:
        kinds = {x.name: x.issym() for x in t.getmembers() if x.isfile() or x.issym()}
    assert kinds == {"results/ecoli/global_exploration/genome_qc/inputs/562.1.fna": True,
                     "results/ecoli/global_exploration/genome_qc/report.csv": False}
    assert row["n_files"] == 2 and row["bytes"] == 2 + len(str(genome.resolve()))
    assert "data/raw/ecoli/genomes" in {r["unit"] for r in rows}


def _upload(tmp_path, env_extra=None):
    fake = tmp_path / "rclone"
    fake.write_text(FAKE_RCLONE)
    fake.chmod(0o755)
    env = {**os.environ, "RCLONE": str(fake), **(env_extra or {})}
    return subprocess.run(["bash", str(PROJECT_ROOT / "scripts" / "backup_upload.sh"),
                           str(tmp_path / "stage"), str(tmp_path / "remote"), "run1",
                           str(tmp_path / "ledger.tsv")], env=env, capture_output=True, text=True)


def test_upload_is_verified_before_the_ledger(tmp_path):
    m = _script()
    _tree(tmp_path)
    rows = m.pack(tmp_path / "stage", tmp_path / "ledger.tsv", current_run=None, threads=1,
                  project=tmp_path)
    archive = rows[0]["archive"]
    r = _upload(tmp_path, {"FAKE_CORRUPT": archive})          # what landed differs
    assert r.returncode == 1 and "do not match" in r.stderr
    assert not (tmp_path / "ledger.tsv").exists() and (tmp_path / "stage" / archive).exists()
    r = _upload(tmp_path, {"FAKE_DUP": archive})              # one name twice on the remote
    assert r.returncode == 1 and "twice" in r.stderr and not (tmp_path / "ledger.tsv").exists()
    r = _upload(tmp_path)
    assert r.returncode == 0, r.stderr
    with open(tmp_path / "ledger.tsv", newline="") as f:
        ledger = list(csv.DictReader(f, delimiter="\t"))
    assert len(ledger) == len(rows) and {x["run"] for x in ledger} == {"run1"}
    assert (tmp_path / "remote" / "ledger.tsv").exists()
    assert not (tmp_path / "stage" / archive).exists()       # staged archives removed
    again = m.pack(tmp_path / "stage", tmp_path / "ledger.tsv", current_run=None, threads=1,
                   project=tmp_path)
    assert again == []
    assert "nothing changed" in _upload(tmp_path).stdout


def test_a_failed_upload_is_tried_again(tmp_path):
    m = _script()
    _tree(tmp_path)
    rows = m.pack(tmp_path / "stage", tmp_path / "ledger.tsv", current_run=None, threads=1,
                  project=tmp_path)
    first, second = rows[0]["archive"], rows[1]["archive"]
    r = _upload(tmp_path, {"FAKE_FAIL_ONCE": first, "BACKUP_ATTEMPTS": "1"})   # no attempt left
    assert r.returncode == 1 and f"FAILED {first}" in r.stderr
    assert not (tmp_path / "ledger.tsv").exists()
    r = _upload(tmp_path, {"FAKE_FAIL_ONCE": second})                          # the second attempt
    assert r.returncode == 0, r.stderr
    assert f"attempt 1 of 3 for {second} ended with exit 1" in r.stderr
    assert (tmp_path / "ledger.tsv").exists()


@pytest.mark.skipif(not shutil.which("timeout"), reason="needs GNU timeout")
def test_a_stalled_upload_is_stopped(tmp_path):
    m = _script()
    _tree(tmp_path)
    rows = m.pack(tmp_path / "stage", tmp_path / "ledger.tsv", current_run=None, threads=1,
                  project=tmp_path)
    archive = rows[0]["archive"]
    r = _upload(tmp_path, {"FAKE_HANG_ONCE": archive, "BACKUP_MIN_SECONDS": "1"})
    assert r.returncode == 0, r.stderr
    assert f"attempt 1 of 3 for {archive} ended with exit 124" in r.stderr

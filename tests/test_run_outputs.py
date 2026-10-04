#!/usr/bin/env python3
"""The checksums of a run's outputs (scripts/run_outputs.py)."""
import hashlib
import importlib.util
import os
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
pytestmark = pytest.mark.unit


def _script():
    spec = importlib.util.spec_from_file_location("amrtest_run_outputs",
                                                  PROJECT_ROOT / "scripts" / "run_outputs.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_outputs_are_located_and_checksums_reused(tmp_path):
    m = _script()
    res = tmp_path / "results"
    config = {"paths_organism": {
        "cv_dir": str(res / "{organism}" / "{antibiotic}" / "cv"),
        "candidates_file": str(res / "{organism}" / "{antibiotic}" / "candidates.csv"),
        "kb_dir": str(res / "kb")}}
    files = {res / "ecoli" / "ampicillin" / "cv" / "oof.csv": b"a,b\n1,2\n",
             res / "ecoli" / "ampicillin" / "candidates.csv": b"pattern_id\n3\n",
             res / "kb" / "kanit.sqlite": b"SQLite",
             res / "notes.txt": b"x"}
    for f, data in files.items():
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_bytes(data)
    run1 = tmp_path / "runs" / "20261005T100000" / "outputs.csv"
    t = m.manifest(config, run1, roots=[str(res)], threads=2).set_index("path")
    oof = str(res / "ecoli" / "ampicillin" / "cv" / "oof.csv")
    assert t.loc[oof, ["location", "organism", "antibiotic"]].tolist() == [
        "cv_dir", "ecoli", "ampicillin"]
    assert t.loc[str(res / "ecoli" / "ampicillin" / "candidates.csv"), "location"] == \
        "candidates_file"
    assert t.loc[str(res / "notes.txt"), "location"] == "other"
    assert t.loc[oof, "sha256"] == hashlib.sha256(b"a,b\n1,2\n").hexdigest()
    run1.parent.mkdir(parents=True)
    t.reset_index().to_csv(run1, index=False)
    # same size and modification time: the earlier checksum is reused, not recomputed
    f = Path(oof)
    st = f.stat()
    f.write_bytes(b"a,b\n9,9\n")
    os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns))
    run2 = tmp_path / "runs" / "20261006T100000" / "outputs.csv"
    t2 = m.manifest(config, run2, roots=[str(res)], threads=2).set_index("path")
    assert t2.loc[oof, "sha256"] == t.loc[oof, "sha256"]
    (res / "kb" / "kanit.sqlite").write_bytes(b"SQLite format 3")
    t3 = m.manifest(config, run2, roots=[str(res)], threads=2).set_index("path")
    assert t3.loc[str(res / "kb" / "kanit.sqlite"), "sha256"] == hashlib.sha256(
        b"SQLite format 3").hexdigest()

#!/usr/bin/env python3
"""03u_unitig_matrix.py: the store's genome set and the unitig-caller call.

The store and model-matrix formats are tested in test_matrix_store.py. Here a
stand-in for unitig-caller (a small script that checks its arguments and writes
an Rtab) exercises the command, its log and the reuse of a finished call.
"""
import os
import stat
import sys

import pandas as pd
import pytest

pytestmark = pytest.mark.unit


@pytest.fixture
def mod(load_script):
    return load_script("03u_unitig_matrix.py")


def test_store_genome_ids_is_the_union_of_the_pairs(mod):
    phenotypes = pd.DataFrame({"Genome ID": ["g1", "g2", "g3", "g4"],
                               "ampicillin": [1, None, 0, 1],
                               "gentamicin": [None, 1, 0, None]})
    qc = pd.DataFrame({"genome_id": ["g1", "g2", "g3", "g4"],
                       "pass_overall": [True, True, True, False]})
    clusters = pd.DataFrame({"Genome ID": ["g1", "g2", "g3", "g4"], "Cluster": [1, 1, 2, 2]})
    # g4 fails QC; g2 is tested only for gentamicin
    assert mod.store_genome_ids(phenotypes, qc, clusters, ["ampicillin"]) == ["g1", "g3"]
    assert mod.store_genome_ids(phenotypes, qc, clusters,
                                ["ampicillin", "gentamicin"]) == ["g1", "g2", "g3"]


def _fake_unitig_caller(tmp_path):
    """Checks the arguments 03u passes and writes <out>.rtab for the refs it got."""
    tool = tmp_path / "unitig-caller"
    tool.write_text(f"""#!{sys.executable}
import sys
from pathlib import Path
a = sys.argv[1:]
assert a[0] == "--call" and "--rtab" in a and a[a.index("--kmer") + 1] == "31"
refs = Path(a[a.index("--refs") + 1]).read_text().split()
names = [Path(r).stem for r in refs]
out = a[a.index("--out") + 1]
print("calling", len(names), "genomes")
with open(out + ".rtab", "w") as f:
    f.write("Unitig_sequence\\t" + "\\t".join(names) + "\\n")
    f.write("ACGT\\t" + "\\t".join("1" for _ in names) + "\\n")
""")
    tool.chmod(tool.stat().st_mode | stat.S_IEXEC)
    return tool


def test_call_unitigs_runs_logs_and_reuses_a_finished_call(mod, tmp_path):
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    for g in ("g1", "g2"):
        (genomes_dir / f"{g}.fna").write_text(">c\nACGT\n")
    tool = _fake_unitig_caller(tmp_path)
    work = tmp_path / "call"
    rtab = mod.call_unitigs(["g1", "g2"], genomes_dir, work, threads=2, k=31, tool=tool)
    assert rtab.read_text().splitlines()[0] == "Unitig_sequence\tg1\tg2"
    assert (work / "unitigs.rtab.done").exists()
    assert "calling 2 genomes" in (work / "unitig_caller.log").read_text()
    refs = (work / "unitig_refs.txt").read_text().split()
    assert refs == [str((genomes_dir / f"{g}.fna").resolve()) for g in ("g1", "g2")]
    # a finished call is reused: the tool is not needed again
    os.remove(tool)
    assert mod.call_unitigs(["g1", "g2"], genomes_dir, work, threads=2, k=31, tool=tool) == rtab


def test_call_unitigs_never_reuses_an_unfinished_rtab(mod, tmp_path):
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    (genomes_dir / "g1.fna").write_text(">c\nACGT\n")
    work = tmp_path / "call"
    work.mkdir()
    (work / "unitigs.rtab").write_text("truncated")       # no .done marker
    tool = _fake_unitig_caller(tmp_path)
    rtab = mod.call_unitigs(["g1"], genomes_dir, work, threads=1, k=31, tool=tool)
    assert rtab.read_text().startswith("Unitig_sequence\tg1")


def test_call_unitigs_stops_on_a_missing_assembly(mod, tmp_path):
    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir()
    with pytest.raises(SystemExit, match="missing"):
        mod.call_unitigs(["g1"], genomes_dir, tmp_path / "call", threads=1, k=31,
                         tool=tmp_path / "unitig-caller")

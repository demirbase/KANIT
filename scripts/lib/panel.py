"""The panel rule: which organism–antibiotic pairs are modelled.

A pair enters the panel when, among the genomes that pass quality control (02d)
and have a lineage (02c), the smaller of the resistant and susceptible classes
holds at least ``min_minority`` genomes. A label that is not a single drug never
enters. A drug that passes the rule but has no class in the antibiotic registry
blocks the panel until the registry is completed, so no pair is decided on an
unclassified drug.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from lib import registry
from lib.config import resolve_path

INCLUDED, EXCLUDED, BLOCKED = "included", "excluded", "blocked"

COLUMNS = ["organism", "antibiotic", "drug_class", "n_tested", "n_eligible",
           "n_resistant", "n_susceptible", "minority", "decision", "reason"]


def input_paths(organism: str, config: dict) -> dict[str, Path]:
    """The three inputs of the panel rule: 00's phenotype matrix, 02d's QC table
    and 02c's lineage table."""
    qc_dir = resolve_path("genome_qc_dir", organism=organism, config=config)
    lineage_dir = resolve_path("lineage_dir", organism=organism, config=config)
    return {
        "phenotypes": resolve_path("metadata_file", organism=organism, config=config),
        "qc_table": qc_dir / f"02d_genome_qc_{organism}.csv",
        "clusters": lineage_dir / "poppunk_clusters.csv",
    }


def read_inputs(organism: str, config: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """(phenotypes, qc_table, clusters) with genome identifiers as strings."""
    paths = input_paths(organism, config)
    return (pd.read_csv(paths["phenotypes"], dtype={"Genome ID": str}, encoding="utf-8"),
            pd.read_csv(paths["qc_table"], dtype={"genome_id": str}, encoding="utf-8"),
            pd.read_csv(paths["clusters"], dtype={"Genome ID": str}, encoding="utf-8"))


def pair_genomes(phenotypes: pd.DataFrame, qc_table: pd.DataFrame, clusters: pd.DataFrame,
                 antibiotic: str) -> pd.DataFrame:
    """Genome ID, label and lineage of the genomes of one pair: eligible and
    tested for the antibiotic, sorted by Genome ID."""
    eligible = eligible_genomes(qc_table, clusters)
    sel = phenotypes.loc[phenotypes["Genome ID"].astype(str).isin(eligible)
                         & phenotypes[antibiotic].notna(), ["Genome ID", antibiotic]]
    lineage = dict(zip(clusters["Genome ID"].astype(str), clusters["Cluster"], strict=True))
    out = pd.DataFrame({"Genome ID": sel["Genome ID"].astype(str),
                        "label": sel[antibiotic].astype(int)})
    out["lineage"] = out["Genome ID"].map(lineage)
    return out.sort_values("Genome ID").reset_index(drop=True)


def decide(n_resistant: int, n_susceptible: int, *, single_drug: bool,
           registered: bool, min_minority: int) -> tuple[str, str]:
    """Decision and reason for one pair. The rules apply in this order."""
    minority = min(n_resistant, n_susceptible)
    if not single_drug:
        return EXCLUDED, "not a single drug"
    if minority < min_minority:
        return EXCLUDED, f"minority {minority} < {min_minority}"
    if not registered:
        return BLOCKED, "passes the rule but has no class in the antibiotic registry"
    return INCLUDED, f"minority {minority} >= {min_minority}"


def eligible_genomes(qc_table: pd.DataFrame, clusters: pd.DataFrame) -> set[str]:
    """Genomes that pass quality control and have a lineage.

    ``qc_table`` is 02d's per-genome table (``genome_id``, ``pass_overall``);
    ``clusters`` is 02c's lineage table (``Genome ID``, ``Cluster``). A blank
    ``pass_overall`` counts as a fail.
    """
    passed = qc_table["pass_overall"].astype(str).str.strip().str.lower().eq("true")
    qc_pass = set(qc_table.loc[passed, "genome_id"].astype(str))
    return qc_pass & set(clusters["Genome ID"].astype(str))


def non_canonical_columns(phenotypes: pd.DataFrame) -> list[tuple[str, str]]:
    """Antibiotic columns whose name is not the registry's canonical spelling, as
    (name, canonical). 00 writes canonical names; any found here mean the matrix
    predates the registry and 00 has to be run again, because merging two
    spellings of one drug needs 00's conflict resolution."""
    out = []
    for ab in phenotypes.columns.drop("Genome ID"):
        canonical = registry.normalize_antibiotic(ab)
        if canonical != ab:
            out.append((ab, canonical))
    return out


def pair_rows(organism: str, phenotypes: pd.DataFrame, eligible: set[str], *,
              min_minority: int) -> list[dict]:
    """One decision row per antibiotic column of a phenotype matrix.

    ``phenotypes`` is 00's matrix: ``Genome ID`` plus one column per antibiotic
    with 1 (resistant), 0 (susceptible) or blank (not tested).
    """
    keep = phenotypes["Genome ID"].astype(str).isin(eligible)
    rows = []
    for ab in phenotypes.columns.drop("Genome ID"):
        col = phenotypes[ab]
        tested = col.notna()
        values = col[keep & tested]
        unexpected = set(values.unique()) - {0, 1}
        if unexpected:
            raise ValueError(f"{organism}/{ab}: phenotype values must be 0 or 1, "
                             f"found {sorted(unexpected)}")
        n_r, n_s = int((values == 1).sum()), int((values == 0).sum())
        drug_class = registry.antibiotic_to_class(ab)
        decision, reason = decide(n_r, n_s, single_drug=registry.is_single_drug(ab),
                                  registered=drug_class is not None,
                                  min_minority=min_minority)
        rows.append({"organism": organism, "antibiotic": ab, "drug_class": drug_class,
                     "n_tested": int(tested.sum()), "n_eligible": int(len(values)),
                     "n_resistant": n_r, "n_susceptible": n_s,
                     "minority": min(n_r, n_s), "decision": decision, "reason": reason})
    return rows


def read_panel(path: str | Path) -> pd.DataFrame:
    """The panel decisions written by 02e."""
    return pd.read_csv(path, encoding="utf-8")


def included_pairs(path: str | Path) -> list[tuple[str, str]]:
    """The (organism, antibiotic) pairs that entered the panel."""
    df = read_panel(path)
    df = df[df["decision"] == INCLUDED]
    return list(zip(df["organism"], df["antibiotic"], strict=True))

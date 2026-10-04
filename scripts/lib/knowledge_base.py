"""KANIT knowledge base v1: schema, identifiers and checks (docs/V1_BILGI_TABANI.md).

The knowledge base stores what the steps measured and inferred and computes
nothing new. Measurements (the evidence layers), inferences (the grades) and
context (NCBI) live in separate tables; no table links context to a grade. It is
built once from the steps' outputs inside one transaction (build_kb.py), checked
by ``validate`` and read-only afterwards.
"""
from __future__ import annotations

import hashlib
import sqlite3

import pandas as pd

from lib import grading
from lib.card_layer import revcomp

SCHEMA_VERSION = "1.0"

SCHEMA = """
-- Release and provenance --------------------------------------------------
CREATE TABLE release (
    kb_version        TEXT NOT NULL,
    schema_version    TEXT NOT NULL,
    created_at        TEXT NOT NULL,
    protocol_version  TEXT NOT NULL,
    protocol_sha256   TEXT NOT NULL,      -- checksum of the protocol text in force
    config_sha256     TEXT NOT NULL,
    code_commit       TEXT,
    code_dirty        INTEGER CHECK (code_dirty IN (0, 1)),
    card_version      TEXT,
    tools             TEXT NOT NULL,      -- JSON {tool: version}
    license           TEXT,               -- data licence (decided in F8.2)
    doi               TEXT
);
CREATE TABLE source_file (              -- every step output the build read
    path    TEXT PRIMARY KEY,
    sha256  TEXT NOT NULL,
    bytes   INTEGER NOT NULL
);
CREATE TABLE data_snapshot (            -- BV-BRC snapshot of each organism (00a, §2.1)
    organism_id   TEXT PRIMARY KEY REFERENCES organism,
    source        TEXT NOT NULL,            -- the API queried
    api_version   TEXT NOT NULL,
    queried_at    TEXT NOT NULL,
    frozen_at     TEXT NOT NULL,
    filters       TEXT NOT NULL,            -- JSON: the genome and the record queries
    n_genomes     INTEGER NOT NULL,         -- phenotyped genomes whose assembly passed
    sha256        TEXT NOT NULL             -- of snapshot.json (file checksums inside)
);
CREATE TABLE reference_database (       -- databases.py manifest (§13)
    name           TEXT PRIMARY KEY,
    version        TEXT NOT NULL,           -- read from the database itself
    downloaded_on  TEXT NOT NULL,
    source         TEXT,
    n_files        INTEGER NOT NULL,
    sha256         TEXT NOT NULL            -- of the sorted per-file checksums
);
CREATE TABLE parameter (                -- every threshold of the protocol
    name              TEXT PRIMARY KEY,
    value             TEXT NOT NULL,
    protocol_section  TEXT NOT NULL
);

-- Reference data ------------------------------------------------------------
CREATE TABLE organism (
    organism_id  TEXT PRIMARY KEY,
    name         TEXT NOT NULL,
    ncbi_taxids  TEXT NOT NULL,      -- its genomes lie at or below these taxa, ';'-joined
    gram_stain   TEXT,
    phylum       TEXT
);
CREATE TABLE antibiotic (
    antibiotic_id      TEXT PRIMARY KEY,  -- canonical registry name
    drug_class         TEXT,              -- registry class; NULL when not registered
    card_drug_classes  TEXT NOT NULL,     -- CARD terms of the class (Appendix A), ';'-joined
    who_aware          TEXT CHECK (who_aware IN ('Access', 'Watch', 'Reserve'))
);
CREATE TABLE aro (                      -- CARD terms hit by RGI (version in release)
    aro_accession  TEXT PRIMARY KEY,
    name           TEXT NOT NULL,
    model_type     TEXT NOT NULL CHECK (model_type IN ('homolog', 'variant', 'overexpression',
                                                       'rrna_variant')),
    gene_family    TEXT,
    drug_class     TEXT,
    mechanism      TEXT
);
CREATE TABLE genome (
    genome_id              TEXT PRIMARY KEY,   -- BV-BRC genome id
    organism_id            TEXT NOT NULL REFERENCES organism,
    checkm2_completeness   REAL,
    checkm2_contamination  REAL,
    n50                    INTEGER,
    n_contigs              INTEGER,
    total_length           INTEGER,
    qc_pass                INTEGER NOT NULL CHECK (qc_pass IN (0, 1)),
    lineage_cluster        TEXT,               -- PopPUNK; NULL when not assigned
    ncbi_taxid             INTEGER NOT NULL,   -- the genome's own taxon
    assembly_accession     TEXT,               -- NCBI identifiers where BV-BRC has them
    sra_accession          TEXT,               -- runs ','-joined
    biosample_accession    TEXT
);
CREATE TABLE phenotype (
    genome_id      TEXT NOT NULL REFERENCES genome,
    antibiotic_id  TEXT NOT NULL REFERENCES antibiotic,
    resistant      INTEGER NOT NULL CHECK (resistant IN (0, 1)),
    PRIMARY KEY (genome_id, antibiotic_id)
);

-- Panel and models ----------------------------------------------------------
CREATE TABLE panel_decision (           -- every organism x antibiotic pair (§3)
    organism_id    TEXT NOT NULL REFERENCES organism,
    antibiotic_id  TEXT NOT NULL,
    decision       TEXT NOT NULL CHECK (decision IN ('included', 'excluded', 'blocked')),
    reason         TEXT,
    n_resistant    INTEGER NOT NULL,
    n_susceptible  INTEGER NOT NULL,
    PRIMARY KEY (organism_id, antibiotic_id)
);
CREATE TABLE model (                    -- one included pair (§3–§6)
    model_id               TEXT PRIMARY KEY,   -- organism__antibiotic
    organism_id            TEXT NOT NULL REFERENCES organism,
    antibiotic_id          TEXT NOT NULL REFERENCES antibiotic,
    n_genomes              INTEGER NOT NULL,
    n_resistant            INTEGER NOT NULL,
    n_susceptible          INTEGER NOT NULL,
    n_lineages             INTEGER NOT NULL,
    largest_lineage_share  REAL NOT NULL,
    n_unitigs              INTEGER NOT NULL,   -- after the frequency filter (§4)
    n_patterns             INTEGER NOT NULL,
    evaluable              INTEGER NOT NULL CHECK (evaluable IN (0, 1)),
    final_n_trees          INTEGER,
    final_params           TEXT,               -- JSON
    UNIQUE (organism_id, antibiotic_id)
);
CREATE TABLE model_genome (
    model_id         TEXT NOT NULL REFERENCES model,
    genome_id        TEXT NOT NULL REFERENCES genome,
    row_index        INTEGER NOT NULL,         -- bit position in pattern.carriers
    resistant        INTEGER NOT NULL CHECK (resistant IN (0, 1)),
    lineage_cluster  TEXT NOT NULL,
    PRIMARY KEY (model_id, genome_id),
    UNIQUE (model_id, row_index)
);
CREATE TABLE model_metric (             -- pooled out-of-fold metrics (§6.3)
    model_id  TEXT NOT NULL REFERENCES model,
    arm       TEXT NOT NULL CHECK (arm IN ('lineage_aware', 'lineage_blind',
                                           'lineage_blind_minus_lineage_aware')),
    metric    TEXT NOT NULL,
    value     REAL,
    ci_low    REAL,                         -- lineage-cluster bootstrap (ROC-AUC only)
    ci_high   REAL,
    PRIMARY KEY (model_id, arm, metric)
);
CREATE TABLE label_permutation (        -- §8.6
    model_id        TEXT PRIMARY KEY REFERENCES model,
    auc_observed    REAL NOT NULL,
    null_mean       REAL NOT NULL,
    null_sd         REAL,
    z               REAL,
    n_permutations  INTEGER NOT NULL,
    p               REAL NOT NULL,
    q               REAL NOT NULL,          -- Benjamini–Hochberg across models
    flag            TEXT NOT NULL CHECK (flag IN ('', 'permutation_not_significant'))
);

-- Features --------------------------------------------------------------------
CREATE TABLE unitig (
    unitig_id  TEXT PRIMARY KEY,            -- KU + SHA-256 of the canonical sequence
    sequence   TEXT NOT NULL UNIQUE,        -- canonical orientation
    length     INTEGER NOT NULL
);
CREATE TABLE pattern (                  -- presence pattern of a model (§4)
    model_id    TEXT NOT NULL REFERENCES model,
    pattern_id  INTEGER NOT NULL,           -- valid within the model only
    n_members   INTEGER NOT NULL,
    n_present   INTEGER NOT NULL,
    carriers    BLOB,                       -- candidates: packed bits in model_genome.row_index order
    PRIMARY KEY (model_id, pattern_id)
);
CREATE TABLE pattern_member (           -- unitigs of the candidate patterns
    model_id    TEXT NOT NULL,
    pattern_id  INTEGER NOT NULL,
    unitig_id   TEXT NOT NULL REFERENCES unitig,
    PRIMARY KEY (model_id, pattern_id, unitig_id),
    FOREIGN KEY (model_id, pattern_id) REFERENCES pattern
);
CREATE TABLE candidate (                -- §7
    model_id    TEXT NOT NULL,
    pattern_id  INTEGER NOT NULL,
    source      TEXT NOT NULL CHECK (source IN ('gain', 'cpss', 'both')),
    gain_rank   INTEGER,
    total_gain  REAL,
    PRIMARY KEY (model_id, pattern_id),
    FOREIGN KEY (model_id, pattern_id) REFERENCES pattern
);

-- Evidence layers (measurements) ---------------------------------------------
CREATE TABLE prevalence_result (        -- §8.2, candidates
    model_id             TEXT NOT NULL,
    pattern_id           INTEGER NOT NULL,
    present_resistant    INTEGER NOT NULL,
    present_susceptible  INTEGER NOT NULL,
    prev_resistant       REAL NOT NULL,
    prev_susceptible     REAL NOT NULL,
    delta                REAL NOT NULL,
    direction            TEXT NOT NULL CHECK (direction IN ('R', 'S', 'none')),
    fisher_p             REAL NOT NULL,
    q                    REAL NOT NULL,
    passes               INTEGER NOT NULL CHECK (passes IN (0, 1)),
    PRIMARY KEY (model_id, pattern_id),
    FOREIGN KEY (model_id, pattern_id) REFERENCES candidate
);
CREATE TABLE mda_result (               -- §8.3, candidates
    model_id                TEXT NOT NULL,
    pattern_id              INTEGER NOT NULL,
    n_fold_models_using     INTEGER NOT NULL,
    auc_observed            REAL NOT NULL,
    auc_permuted            REAL NOT NULL,
    mda                     REAL NOT NULL,
    n_permuted_ge_observed  INTEGER NOT NULL,
    p                       REAL NOT NULL,
    q                       REAL NOT NULL,
    passes                  INTEGER NOT NULL CHECK (passes IN (0, 1)),
    PRIMARY KEY (model_id, pattern_id),
    FOREIGN KEY (model_id, pattern_id) REFERENCES candidate
);
CREATE TABLE mda_cluster (              -- §8.3 sensitivity analysis, candidates
    model_id          TEXT NOT NULL,
    pattern_id        INTEGER NOT NULL,
    cluster           INTEGER NOT NULL,
    cluster_size      INTEGER NOT NULL,
    cluster_patterns  TEXT NOT NULL,
    mda               REAL NOT NULL,
    p                 REAL NOT NULL,
    q                 REAL NOT NULL,
    passes            INTEGER NOT NULL CHECK (passes IN (0, 1)),
    PRIMARY KEY (model_id, pattern_id),
    FOREIGN KEY (model_id, pattern_id) REFERENCES candidate
);
CREATE TABLE cpss_result (              -- §8.4, prefilter and candidates
    model_id      TEXT NOT NULL,
    pattern_id    INTEGER NOT NULL,
    in_prefilter  INTEGER NOT NULL CHECK (in_prefilter IN (0, 1)),
    chi2          REAL,
    n_selected    INTEGER NOT NULL,
    pi            REAL NOT NULL,
    passes        INTEGER NOT NULL CHECK (passes IN (0, 1)),
    PRIMARY KEY (model_id, pattern_id),
    FOREIGN KEY (model_id, pattern_id) REFERENCES pattern
);
CREATE TABLE pyseer_result (            -- §8.5, every tested pattern
    model_id    TEXT NOT NULL,
    pattern_id  INTEGER NOT NULL,
    af          REAL,
    beta        REAL,
    beta_se     REAL,
    lrt_p       REAL,
    notes       TEXT NOT NULL,
    passes      INTEGER NOT NULL CHECK (passes IN (0, 1)),
    PRIMARY KEY (model_id, pattern_id),
    FOREIGN KEY (model_id, pattern_id) REFERENCES pattern
);
CREATE TABLE card_hit (                 -- §8.1: RGI hits overlapping a candidate unitig
    model_id       TEXT NOT NULL REFERENCES model,
    unitig_id      TEXT NOT NULL REFERENCES unitig,
    aro_accession  TEXT NOT NULL REFERENCES aro,
    PRIMARY KEY (model_id, unitig_id, aro_accession)
);
CREATE TABLE card_annotation (          -- §8.1: CARD state of every candidate unitig
    model_id         TEXT NOT NULL REFERENCES model,
    unitig_id        TEXT NOT NULL REFERENCES unitig,
    mode             TEXT NOT NULL CHECK (mode IN ('allele_aware', 'homolog_only')),
    state            TEXT NOT NULL CHECK (state IN ('b', 'card_hit_without_b', 'no_card_hit')),
    reasons          TEXT NOT NULL,
    n_b              INTEGER NOT NULL,
    located_genomes  TEXT NOT NULL,
    PRIMARY KEY (model_id, unitig_id, mode)
);

-- Inference ---------------------------------------------------------------------
CREATE TABLE grade (                    -- §9, both rules
    model_id       TEXT NOT NULL,
    pattern_id     INTEGER NOT NULL,
    rule           TEXT NOT NULL CHECK (rule IN ('allele_aware', 'homolog_only')),
    primary_rule   INTEGER NOT NULL CHECK (primary_rule IN (0, 1)),
    card_state     TEXT NOT NULL CHECK (card_state IN ('b', 'card_hit_without_b', 'no_card_hit')),
    card_reasons   TEXT NOT NULL,
    n_layers       INTEGER NOT NULL CHECK (n_layers BETWEEN 0 AND 4),
    layers_passed  TEXT NOT NULL,
    grade          TEXT NOT NULL CHECK (grade IN ('confirmed', 'strong_novel', 'candidate',
                                                  'weak', 'none')),
    PRIMARY KEY (model_id, pattern_id, rule),
    FOREIGN KEY (model_id, pattern_id) REFERENCES candidate
);

-- Context (not evidence; §12) -------------------------------------------------------
CREATE TABLE unitig_context (           -- NCBI nt, searched within the organism
    unitig_id       TEXT NOT NULL REFERENCES unitig,
    organism_id     TEXT NOT NULL REFERENCES organism,
    source          TEXT NOT NULL,          -- ncbi_nt_remote
    queried_on      TEXT NOT NULL,
    nt_release      TEXT,                   -- 'Posted date' of the BLAST report
    n_hits          INTEGER NOT NULL,
    best_accession  TEXT,
    best_title      TEXT,
    best_identity   REAL,
    best_coverage   REAL,
    best_evalue     REAL,
    gene            TEXT,
    product         TEXT,
    plasmid_share   REAL,                   -- share of the hits on plasmid records
    PRIMARY KEY (unitig_id, organism_id, source)
);

-- Comparison with genotype-based prediction (§11) ------------------------------------
CREATE TABLE external_call (            -- AMRFinderPlus calls, the reference of H7 and H6
    genome_id       TEXT NOT NULL REFERENCES genome,
    call_index      INTEGER NOT NULL,
    element_symbol  TEXT NOT NULL,
    type            TEXT NOT NULL,
    subtype         TEXT NOT NULL,
    scope           TEXT NOT NULL,
    class           TEXT NOT NULL,
    subclass        TEXT NOT NULL,
    PRIMARY KEY (genome_id, call_index)
);
CREATE TABLE external_comparison (
    model_id               TEXT NOT NULL REFERENCES model,
    tool                   TEXT NOT NULL CHECK (tool IN ('amrfinderplus', 'resfinder', 'rgi_all',
                                                         'rgi_without_near_universal', 'model')),
    assessable             INTEGER NOT NULL CHECK (assessable IN (0, 1)),
    n                      INTEGER,
    n_resistant            INTEGER,
    tp                     INTEGER,
    fp                     INTEGER,
    tn                     INTEGER,
    fn                     INTEGER,
    sensitivity            REAL,
    specificity            REAL,
    balanced_accuracy      REAL,
    very_major_error_rate  REAL,
    major_error_rate       REAL,
    PRIMARY KEY (model_id, tool)
);

-- Views -----------------------------------------------------------------------------
CREATE VIEW v_biomarker AS
SELECT c.model_id, m.organism_id, m.antibiotic_id, c.pattern_id, p.n_members, p.n_present,
       c.source AS candidate_source, c.gain_rank,
       gp.rule AS grading_rule, gp.grade,
       ga.grade AS grade_allele_aware, gh.grade AS grade_homolog_only,
       ga.card_state AS card_state_allele_aware, ga.card_reasons AS card_reasons_allele_aware,
       gh.card_state AS card_state_homolog_only,
       ga.n_layers, ga.layers_passed,
       pr.delta, pr.direction, pr.q AS prevalence_q, md.mda, md.q AS mda_q, cs.pi,
       py.lrt_p AS pyseer_p,
       (SELECT group_concat(DISTINCT a.name)
          FROM pattern_member pm
          JOIN card_hit h ON h.model_id = pm.model_id AND h.unitig_id = pm.unitig_id
          JOIN aro a ON a.aro_accession = h.aro_accession
         WHERE pm.model_id = c.model_id AND pm.pattern_id = c.pattern_id) AS card_genes
  FROM candidate c
  JOIN model m ON m.model_id = c.model_id
  JOIN pattern p ON p.model_id = c.model_id AND p.pattern_id = c.pattern_id
  JOIN grade gp ON gp.model_id = c.model_id AND gp.pattern_id = c.pattern_id AND gp.primary_rule = 1
  JOIN grade ga ON ga.model_id = c.model_id AND ga.pattern_id = c.pattern_id AND ga.rule = 'allele_aware'
  JOIN grade gh ON gh.model_id = c.model_id AND gh.pattern_id = c.pattern_id AND gh.rule = 'homolog_only'
  JOIN prevalence_result pr ON pr.model_id = c.model_id AND pr.pattern_id = c.pattern_id
  JOIN mda_result md ON md.model_id = c.model_id AND md.pattern_id = c.pattern_id
  JOIN cpss_result cs ON cs.model_id = c.model_id AND cs.pattern_id = c.pattern_id
  JOIN pyseer_result py ON py.model_id = c.model_id AND py.pattern_id = c.pattern_id;

CREATE VIEW v_unitig AS
SELECT u.unitig_id, u.sequence, u.length, pm.model_id, pm.pattern_id, g.grade,
       ca.state AS card_state_allele_aware, ca.reasons AS card_reasons_allele_aware
  FROM pattern_member pm
  JOIN unitig u ON u.unitig_id = pm.unitig_id
  JOIN grade g ON g.model_id = pm.model_id AND g.pattern_id = pm.pattern_id AND g.primary_rule = 1
  JOIN card_annotation ca ON ca.model_id = pm.model_id AND ca.unitig_id = pm.unitig_id
                         AND ca.mode = 'allele_aware';

CREATE VIEW v_model AS
SELECT m.model_id, m.organism_id, m.antibiotic_id, m.n_genomes, m.n_resistant, m.n_susceptible,
       m.n_lineages, m.evaluable,
       aw.value AS roc_auc_lineage_aware, aw.ci_low AS roc_auc_lineage_aware_low,
       aw.ci_high AS roc_auc_lineage_aware_high, bl.value AS roc_auc_lineage_blind,
       lp.z AS label_permutation_z, lp.q AS label_permutation_q, lp.flag,
       (SELECT count(*) FROM candidate c WHERE c.model_id = m.model_id) AS n_candidates,
       (SELECT count(*) FROM grade g WHERE g.model_id = m.model_id AND g.primary_rule = 1
                                       AND g.grade = 'confirmed') AS n_confirmed,
       (SELECT count(*) FROM grade g WHERE g.model_id = m.model_id AND g.primary_rule = 1
                                       AND g.grade = 'strong_novel') AS n_strong_novel
  FROM model m
  LEFT JOIN model_metric aw ON aw.model_id = m.model_id AND aw.arm = 'lineage_aware'
                           AND aw.metric = 'roc_auc'
  LEFT JOIN model_metric bl ON bl.model_id = m.model_id AND bl.arm = 'lineage_blind'
                           AND bl.metric = 'roc_auc'
  LEFT JOIN label_permutation lp ON lp.model_id = m.model_id;

CREATE VIEW v_novel AS                  -- context is shown beside, never as evidence
SELECT n.*, c.best_title AS context_best_title, c.gene AS context_gene,
       c.product AS context_product, c.plasmid_share AS context_plasmid_share,
       0 AS context_is_evidence
  FROM (SELECT b.model_id, b.organism_id, b.pattern_id, b.n_members, b.layers_passed,
               (SELECT pm.unitig_id FROM pattern_member pm JOIN unitig u ON u.unitig_id = pm.unitig_id
                 WHERE pm.model_id = b.model_id AND pm.pattern_id = b.pattern_id
                 ORDER BY u.length DESC, u.unitig_id LIMIT 1) AS longest_unitig_id
          FROM v_biomarker b
         WHERE b.grade = 'strong_novel') n
  LEFT JOIN unitig_context c ON c.unitig_id = n.longest_unitig_id
                            AND c.organism_id = n.organism_id;
"""

# The protocol's thresholds: (parameter name, config path, protocol section).
PARAMETERS = (
    ("panel.min_minority", ("panel", "min_minority"), "§3"),
    ("unitig.k", ("unitig", "k"), "§4"),
    ("unitig.min_support", ("unitig", "min_support"), "§4"),
    ("hpo.n_trials", ("hpo", "n_trials"), "§5"),
    ("hpo.search_max_genomes", ("hpo", "search_max_genomes"), "§5"),
    ("hpo.early_stopping_rounds", ("hpo", "early_stopping_rounds"), "§5"),
    ("hpo.max_rounds", ("hpo", "max_rounds"), "§5"),
    ("cv.n_repeats", ("cv", "n_repeats"), "§6.1"),
    ("cv.n_folds", ("cv", "n_folds"), "§6.1"),
    ("cv.min_minority_per_test_fold", ("cv", "min_minority_per_test_fold"), "§6.2"),
    ("cv.threshold", ("cv", "threshold"), "§6.3"),
    ("cv.n_bootstrap", ("cv", "n_bootstrap"), "§6.3"),
    ("candidates.top_gain", ("candidates", "top_gain"), "§7"),
    ("card.near_universal", ("card", "near_universal"), "§8.1"),
    ("card.max_located_genomes", ("card", "max_located_genomes"), "§8.1"),
    ("card.min_overlap", ("card", "min_overlap"), "§8.1"),
    ("prevalence.min_delta", ("prevalence", "min_delta"), "§8.2"),
    ("prevalence.alpha", ("prevalence", "alpha"), "§8.2"),
    ("mda.n_permutations", ("mda", "n_permutations"), "§8.3"),
    ("mda.alpha", ("mda", "alpha"), "§8.3"),
    ("mda.cluster_r", ("mda", "cluster_r"), "§8.3"),
    ("cpss.prefilter", ("cpss", "prefilter"), "§8.4"),
    ("cpss.n_pairs", ("cpss", "n_pairs"), "§8.4"),
    ("cpss.q", ("cpss", "q"), "§8.4"),
    ("cpss.pi_threshold", ("cpss", "pi_threshold"), "§8.4"),
    ("pyseer.kinship_every", ("pyseer", "kinship_every"), "§8.5"),
    ("pyseer.alpha", ("pyseer", "alpha"), "§8.5"),
    ("label_permutation.n_permutations", ("label_permutation", "n_permutations"), "§8.6"),
    ("label_permutation.alpha", ("label_permutation", "alpha"), "§8.6"),
    ("grading.rule", ("grading", "rule"), "§9"),
)


def unitig_id(sequence: str) -> tuple[str, str]:
    """(identifier, canonical sequence): the canonical orientation is the smaller of
    the sequence and its reverse complement; the identifier is KU and the first 20
    hexadecimal digits of its SHA-256."""
    s = sequence.strip().upper()
    canonical = min(s, revcomp(s))
    return "KU" + hashlib.sha256(canonical.encode("ascii")).hexdigest()[:20], canonical


def parameters(config: dict) -> list[tuple[str, str, str]]:
    rows = []
    for name, (section, key), where in PARAMETERS:
        rows.append((name, str(config[section][key]), where))
    return rows


def create(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA foreign_keys = ON")
    conn.executescript(SCHEMA)


def insert(conn: sqlite3.Connection, table: str, df: pd.DataFrame) -> int:
    """Insert every row of ``df`` (its columns must be columns of ``table``)."""
    if df.empty:
        return 0
    cols = list(df.columns)
    sql = f"INSERT INTO {table} ({', '.join(cols)}) VALUES ({', '.join('?' * len(cols))})"
    rows = [tuple(None if pd.isna(v) else (v.item() if hasattr(v, "item") else v) for v in r)
            for r in df.itertuples(index=False, name=None)]
    conn.executemany(sql, rows)
    return len(rows)


def validate(conn: sqlite3.Connection) -> dict:
    """Structural checks of a built knowledge base; raises on the first failure."""
    problems = []
    if conn.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
        problems.append("SQLite integrity check failed")
    fk = conn.execute("PRAGMA foreign_key_check").fetchall()
    if fk:
        problems.append(f"{len(fk)} foreign key violation(s), e.g. {fk[:3]}")
    for table, n in (("prevalence_result", 1), ("mda_result", 1), ("cpss_result", 1),
                     ("pyseer_result", 1), ("mda_cluster", 1), ("grade", 2)):
        bad = conn.execute(
            f"SELECT count(*) FROM candidate c WHERE (SELECT count(*) FROM {table} t "
            f"WHERE t.model_id = c.model_id AND t.pattern_id = c.pattern_id) != {n}").fetchone()[0]
        if bad:
            problems.append(f"{bad} candidate(s) without exactly {n} row(s) in {table}")
    bad = conn.execute(
        "SELECT count(*) FROM candidate c WHERE (SELECT count(*) FROM grade g WHERE "
        "g.model_id = c.model_id AND g.pattern_id = c.pattern_id AND g.primary_rule = 1) != 1"
    ).fetchone()[0]
    if bad:
        problems.append(f"{bad} candidate(s) without exactly one primary grade")
    bad = conn.execute("SELECT count(*) FROM organism o WHERE NOT EXISTS (SELECT 1 FROM "
                       "data_snapshot d WHERE d.organism_id = o.organism_id)").fetchone()[0]
    if bad:
        problems.append(f"{bad} organism(s) without a data snapshot")
    bad = conn.execute("SELECT count(*) FROM model m WHERE (SELECT count(*) FROM "
                       "external_comparison e WHERE e.model_id = m.model_id) != 5").fetchone()[0]
    if bad:
        problems.append(f"{bad} model(s) without the five rows of the external comparison")
    missing = conn.execute(
        "SELECT count(*) FROM pattern_member pm WHERE NOT EXISTS (SELECT 1 FROM card_annotation ca "
        "WHERE ca.model_id = pm.model_id AND ca.unitig_id = pm.unitig_id)").fetchone()[0]
    if missing:
        problems.append(f"{missing} candidate unitig(s) without a CARD annotation")
    # the grades follow from the stored layers and CARD states (lib.grading)
    q = """SELECT g.model_id, g.pattern_id, g.rule, g.card_state, g.grade, g.n_layers,
                  pr.passes, md.passes, cs.passes, py.passes
             FROM grade g
             JOIN prevalence_result pr USING (model_id, pattern_id)
             JOIN mda_result md USING (model_id, pattern_id)
             JOIN cpss_result cs USING (model_id, pattern_id)
             JOIN pyseer_result py USING (model_id, pattern_id)"""
    n_checked = 0
    for mid, pid, rule, state, g, n_layers, *passes in conn.execute(q):
        layers = dict(zip(grading.LAYERS, (bool(x) for x in passes), strict=True))
        if grading.grade(state, layers) != g or grading.n_passed(layers) != n_layers:
            problems.append(f"grade of {mid} pattern {pid} ({rule}) does not follow its layers")
            break
        n_checked += 1
    if problems:
        raise ValueError("knowledge base check failed: " + "; ".join(problems))
    counts = {t: conn.execute(f"SELECT count(*) FROM {t}").fetchone()[0]
              for (t,) in conn.execute("SELECT name FROM sqlite_master WHERE type = 'table' "
                                       "ORDER BY name")}
    return {"tables": counts, "grades_rechecked": n_checked}

#!/usr/bin/env nextflow
/*
 * KANIT v1.0 — evidence-graded knowledge base of AMR biomarkers (protocol
 * docs/V1_PROTOKOL.md; every analysis setting in config/config.yaml).
 *
 * Entries, run in this order:
 *   -entry DOWNLOAD   login node (internet): genomes and phenotypes from BV-BRC
 *   (default)         quality control, lineages, panel, features, nested
 *                     cross-validation, evidence layers, grades, label
 *                     permutation and the comparison with genotype-based tools
 *   -entry CONTEXT    login node (internet): NCBI context of the candidate unitigs
 *   -entry KB         the knowledge base and the hypothesis tests
 *
 * The steps write into the shared results tree; every task emits a receipt that
 * the tasks depending on it take as input.
 *
 *   nextflow run main.nf -profile truba --organisms ecoli,kpneumoniae
 *   nextflow run main.nf -stub-run -profile test --organisms ecoli \
 *       --stub_panel tests/data/stub_panel_decisions.csv
 */

include { validateParameters } from 'plugin/nf-schema'
include { GENOMES } from './subworkflows/local/genomes'
include { MODELS } from './subworkflows/local/models'
include { EVIDENCE } from './subworkflows/local/evidence'
include { COMPARISON } from './subworkflows/local/comparison'
include { DOWNLOAD_BVBRC; PREPARE_METADATA } from './modules/local/download'
include { CONTEXT_QUERY; CONTEXT_BUILD } from './modules/local/context'
include { BUILD_KB; HYPOTHESES } from './modules/local/kb'

// the organisms asked for, every one of them in the registry
def organisms() {
    validateParameters()
    def registry = new org.yaml.snakeyaml.Yaml().load(
        file("${projectDir}/config/registry/organisms.yaml").text)
    def known = registry.organisms.keySet()
    def wanted = params.organisms.tokenize(',')
    def unknown = wanted.findAll { !(it in known) }
    if (unknown) {
        error "Unknown organism(s) ${unknown.join(', ')}; the registry has ${known.sort().join(', ')}"
    }
    if (workflow.stubRun && !params.stub_panel) {
        error "-stub-run needs --stub_panel (the panel decisions to fan the models out from)"
    }
    wanted
}

workflow {
    def orgs = Channel.fromList(organisms())
    GENOMES(orgs)
    modelled = GENOMES.out.pairs.map { it.organism }.unique()
    MODELS(GENOMES.out.pairs, GENOMES.out.stores)
    EVIDENCE(GENOMES.out.stores.join(modelled.map { [it] }), MODELS.out.folds, MODELS.out.units,
             MODELS.out.finals)
    COMPARISON(modelled, GENOMES.out.panel, EVIDENCE.out.rgi, MODELS.out.metrics)
}

workflow DOWNLOAD {
    def orgs = Channel.fromList(organisms())
    DOWNLOAD_BVBRC(orgs)
    PREPARE_METADATA(DOWNLOAD_BVBRC.out.done)
}

workflow CONTEXT {
    def orgs = Channel.fromList(organisms())
    CONTEXT_QUERY(orgs)
    CONTEXT_BUILD(CONTEXT_QUERY.out.done.map { org, r -> org })
}

workflow KB {
    def orgs = Channel.fromList(organisms())
    CONTEXT_BUILD(orgs)
    BUILD_KB(CONTEXT_BUILD.out.done.map { org, r -> r }.collect())
    HYPOTHESES(BUILD_KB.out.done)
}

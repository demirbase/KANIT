// Every organism: quality control and lineages, then the panel of all organisms and,
// after it, the unitig store of each organism with panel pairs. Emits the panel pairs.
// An organism's tasks after the panel take its pairs, not the panel's receipt, so that a
// run with more organisms leaves the earlier ones cached.

include { QC_PREP; CHECKM2; QUAST; QC_POST; LINEAGE; PANEL; STORE } from '../../modules/local/genomes'
include { meta } from '../../modules/local/common'

workflow GENOMES {
    take:
    organisms                       // organism ids

    main:
    QC_PREP(organisms)
    CHECKM2(QC_PREP.out.paths)
    QUAST(QC_PREP.out.paths)
    QC_POST(CHECKM2.out.done.join(QUAST.out.done).map { org, a, b -> [org, [a, b]] })
    LINEAGE(organisms)
    ready = QC_POST.out.done.join(LINEAGE.out.done).map { org, a, b -> [org, [a, b]] }
    PANEL(organisms.collect(), ready.flatMap { org, rs -> rs }.collect())
    def wanted = params.organisms.tokenize(',') as Set
    pairs = PANEL.out.decisions
        .splitCsv(header: true)
        .filter { row -> row.decision == 'included' && row.organism in wanted }
        .map { row -> meta(row.organism, row.antibiotic) }
    // every organism with pairs: the antibiotics of its pairs and its QC and lineage
    // receipts. The store (03u) and the genotype-based tools (16) read the panel through
    // them, not through the panel's receipt: the panel is rerun whenever the organisms of a
    // run change, while the pairs of an organism stay the same, so its tasks stay cached
    prepared = pairs.map { [it.organism, it.antibiotic] }.groupTuple()
        .map { org, abs -> [org, abs.sort(false)] }
        .join(ready)
    STORE(prepared)

    emit:
    pairs    = pairs                // meta of every panel pair
    stores   = STORE.out.done       // organism, receipt
    prepared = prepared             // organism, its pairs' antibiotics, QC and lineage receipts
}

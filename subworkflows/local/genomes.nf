// Every organism: quality control and lineages, then the panel of all organisms and,
// after it, the unitig store of each organism with panel pairs. Emits the panel pairs.

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
    // the store holds the genomes of the organism's panel pairs (03u reads the panel):
    // it waits for the panel and is built only for the organisms that have pairs
    withPairs = pairs.map { it.organism }.unique().map { [it] }
    STORE(ready.join(withPairs).combine(PANEL.out.receipt).map { org, rs, p -> [org, rs + [p]] })

    emit:
    pairs  = pairs                  // meta of every panel pair
    stores = STORE.out.done         // organism, receipt
    panel  = PANEL.out.decisions    // the panel decisions
}

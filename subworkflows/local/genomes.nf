// Every organism: quality control and lineages, then the panel of all organisms and
// each organism's unitig store. Emits the panel pairs.

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
    STORE(ready)
    def wanted = params.organisms.tokenize(',') as Set
    pairs = PANEL.out.decisions
        .splitCsv(header: true)
        .filter { row -> row.decision == 'included' && row.organism in wanted }
        .map { row -> meta(row.organism, row.antibiotic) }

    emit:
    pairs  = pairs                  // meta of every panel pair
    stores = STORE.out.done         // organism, receipt
    panel  = PANEL.out.decisions    // the panel decisions
}

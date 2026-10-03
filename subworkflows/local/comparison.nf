// Every organism with panel pairs: AMRFinderPlus and ResFinder on its genomes and the
// comparison of the genotype-based predictions with the models (16).

include { EXTERNAL_PREP; EXTERNAL_RUN; EXTERNAL_COLLECT; EXTERNAL_COMPARE } from '../../modules/local/external'

workflow COMPARISON {
    take:
    organisms                       // organism ids with panel pairs
    panel                           // the panel decisions
    rgi                             // organism, RGI receipt
    metrics                         // meta, metrics receipt

    main:
    EXTERNAL_PREP(organisms.combine(panel))
    EXTERNAL_RUN(EXTERNAL_PREP.out.script)
    EXTERNAL_COLLECT(EXTERNAL_RUN.out.done)
    perOrganism = metrics.map { m, r -> [m.organism, r] }.groupTuple()
    EXTERNAL_COMPARE(EXTERNAL_COLLECT.out.done.join(rgi).join(perOrganism, remainder: true)
                           .map { org, a, b, cs -> [org, [a, b] + (cs ?: [])] })

    emit:
    compare = EXTERNAL_COMPARE.out.done
}

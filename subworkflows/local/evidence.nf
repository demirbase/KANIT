// Every evaluable pair: RGI on its organism, stability selection and the candidates,
// the evidence layers, the CARD layer, the grades and the label permutation test.

include { RGI_LOAD; RGI_RUN; RGI_COLLECT; CPSS_PREFILTER; CPSS_CHUNK; CPSS_SELECT; PREVALENCE;
          MDA; CARD_LAYER; PYSEER_PREP; PYSEER_LMM; PYSEER_POST; GRADING; LP_CHUNK; LP_METRICS;
          LP_ACROSS } from '../../modules/local/evidence'
include { kanitConfig } from '../../modules/local/common'

workflow EVIDENCE {
    take:
    stores                          // organism, store receipt (organisms with pairs)
    folds                           // meta, folds receipt
    units                           // meta, outer fold receipts
    finals                          // meta, final model receipt

    main:
    def cfg = kanitConfig()
    def nShards = params.rgi_shards as int
    def nCpss = Math.ceil(cfg.cpss.n_pairs / cfg.cpss.chunk) as int
    def nLp = Math.ceil(cfg.label_permutation.n_permutations / cfg.label_permutation.chunk) as int
    def nFolds = cfg.cv.n_folds as int

    RGI_LOAD()
    RGI_RUN(stores.combine(RGI_LOAD.out.done)
                  .flatMap { org, s, l -> (0..<nShards).collect { k -> [org, k, nShards, [s, l]] } })
    RGI_COLLECT(RGI_RUN.out.done.groupTuple(size: nShards))

    CPSS_PREFILTER(folds)
    CPSS_CHUNK(CPSS_PREFILTER.out.done.join(finals)
                     .flatMap { m, a, b -> (0..<nCpss).collect { c -> [m, c, [a, b]] } })
    CPSS_SELECT(CPSS_CHUNK.out.done.groupTuple(size: nCpss).join(finals)
                      .map { m, rs, f -> [m, rs + [f]] })
    candidates = CPSS_SELECT.out.done

    PREVALENCE(candidates)
    MDA(candidates.join(units).map { m, c, us -> [m, [c] + us] })
    CARD_LAYER(candidates.map { m, c -> [m.organism, m, c] }.combine(RGI_COLLECT.out.done, by: 0)
                     .map { org, m, c, g -> [m, [c, g]] })
    PYSEER_PREP(candidates)
    PYSEER_LMM(PYSEER_PREP.out.script)
    PYSEER_POST(PYSEER_LMM.out.done)
    GRADING(CARD_LAYER.out.done.join(PREVALENCE.out.done).join(MDA.out.done)
                  .join(PYSEER_POST.out.done).map { m, a, b, c, d -> [m, [a, b, c, d]] })

    LP_CHUNK(units.flatMap { m, us ->
        (0..<nFolds).collectMany { k -> (0..<nLp).collect { c -> [m, k, c, us] } } })
    LP_METRICS(LP_CHUNK.out.done.groupTuple(size: nFolds * nLp))
    LP_ACROSS(LP_METRICS.out.done.map { m, r -> r }.collect())

    emit:
    grading = GRADING.out.done      // meta, grades receipt
    rgi     = RGI_COLLECT.out.done  // organism, RGI receipt
    lp      = LP_ACROSS.out.done    // the across-models table's receipt
}

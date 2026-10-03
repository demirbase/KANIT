// Every panel pair: the model matrix, the outer folds, every evaluable outer fold,
// the final model and the out-of-fold metrics. A pair that is not evaluable stops
// after its folds (protocol §6.2).

include { MODEL_MATRIX; CV_FOLDS; CV_UNIT; CV_FINAL; CV_METRICS } from '../../modules/local/models'

workflow MODELS {
    take:
    pairs                           // meta of every panel pair
    stores                          // organism, store receipt

    main:
    MODEL_MATRIX(pairs.map { m -> [m.organism, m] }.combine(stores, by: 0)
                      .map { org, m, r -> [m, r] })
    CV_FOLDS(MODEL_MATRIX.out.done)
    evaluable = CV_FOLDS.out.units.filter { m, u, r, e -> e == 'true' }
    units = evaluable.flatMap { m, u, r, e ->
        def lines = u.readLines().findAll { it.trim() }
        lines.collect { l ->
            def (arm, repeat, fold) = l.tokenize(' ')
            [m, lines.size(), arm, repeat as int, fold as int, r]
        }
    }
    CV_UNIT(units)
    unitsDone = CV_UNIT.out.done
        .map { m, n, r -> [groupKey(m, n), r] }
        .groupTuple()
        .map { k, rs -> [k.getGroupTarget(), rs] }
    folds = evaluable.map { m, u, r, e -> [m, r] }
    CV_FINAL(folds)
    CV_METRICS(unitsDone)

    emit:
    folds   = folds                 // meta, folds receipt (evaluable pairs)
    units   = unitsDone             // meta, every outer fold's receipt
    finals  = CV_FINAL.out.done     // meta, final model receipt
    metrics = CV_METRICS.out.done   // meta, metrics receipt
}

// One panel pair: its model matrix (03u) and the nested cross-validation with the
// final model (04).

include { py; receipt; stubReceipt; kanitConfig } from './common'

process MODEL_MATRIX {
    tag "${meta.id}"
    label 'matrix'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('03u_unitig_matrix.py')} --organism ${meta.organism} --antibiotic ${meta.antibiotic} \\
        --threads ${task.cpus}
    ${receipt(task, 'model_matrix', meta.id)}
    """

    stub:
    stubReceipt(task, 'model_matrix', meta.id)
}

process CV_FOLDS {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('units.txt'), path('receipt.json'), env(EVALUABLE), emit: units

    script:
    def args = "--organism ${meta.organism} --antibiotic ${meta.antibiotic}"
    """
    ${py('04_nested_cv.py')} folds ${args}
    ${py('04_nested_cv.py')} units ${args} > units.txt
    EVALUABLE=\$( [ -s units.txt ] && echo true || echo false )
    ${receipt(task, 'cv_folds', meta.id)}
    """

    stub:
    def cv = kanitConfig().cv
    def lines = ['lineage_aware', 'lineage_blind'].collectMany { a ->
        (1..cv.n_repeats).collectMany { r -> (0..<cv.n_folds).collect { k -> "${a} ${r} ${k}" } } }
    """
    printf '%s\\n' ${lines.collect { "'${it}'" }.join(' ')} > units.txt
    EVALUABLE=true
    ${stubReceipt(task, 'cv_folds', meta.id)}
    """
}

process CV_UNIT {
    tag "${meta.id} ${arm} r${repeat} f${fold}"
    label 'cv_unit'

    input:
    tuple val(meta), val(n), val(arm), val(repeat), val(fold), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), val(n), path('receipt.json'), emit: done

    script:
    """
    ${py('04_nested_cv.py')} unit --organism ${meta.organism} --antibiotic ${meta.antibiotic} \\
        --arm ${arm} --repeat ${repeat} --fold ${fold} --threads ${task.cpus}
    ${receipt(task, 'cv_unit', "${meta.id} ${arm} ${repeat} ${fold}")}
    """

    stub:
    stubReceipt(task, 'cv_unit', "${meta.id} ${arm} ${repeat} ${fold}")
}

process CV_FINAL {
    tag "${meta.id}"
    label 'cv_unit'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('04_nested_cv.py')} final --organism ${meta.organism} --antibiotic ${meta.antibiotic} \\
        --threads ${task.cpus}
    ${receipt(task, 'cv_final', meta.id)}
    """

    stub:
    stubReceipt(task, 'cv_final', meta.id)
}

process CV_METRICS {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('04_nested_cv.py')} metrics --organism ${meta.organism} --antibiotic ${meta.antibiotic}
    ${receipt(task, 'cv_metrics', meta.id)}
    """

    stub:
    stubReceipt(task, 'cv_metrics', meta.id)
}

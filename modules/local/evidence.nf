// Evidence of one panel pair: RGI on the organism's genomes (08), stability selection
// and the candidates (13), prevalence (10), MDA (12), the CARD layer (09), pyseer
// (14), the grades (14b) and the label permutation test (12b).

include { py; receipt; stubReceipt } from './common'

def pair(meta) {
    "--organism ${meta.organism} --antibiotic ${meta.antibiotic}"
}

process RGI_LOAD {
    label 'rgi'

    output:
    path 'receipt.json', emit: done

    script:
    """
    ${py('08_rgi.py')} load
    ${receipt('rgi_load', 'card')}
    """

    stub:
    stubReceipt('rgi_load', 'card')
}

process RGI_RUN {
    tag "${org} ${shard}/${n}"
    label 'rgi'

    input:
    tuple val(org), val(shard), val(n), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('08_rgi.py')} run --organism ${org} --shard ${shard}/${n} --threads ${task.cpus}
    ${receipt('rgi_run', "${org} ${shard}/${n}")}
    """

    stub:
    stubReceipt('rgi_run', "${org} ${shard}/${n}")
}

process RGI_COLLECT {
    tag "${org}"
    label 'rgi'

    input:
    tuple val(org), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('08_rgi.py')} collect --organism ${org}
    ${receipt('rgi_collect', org)}
    """

    stub:
    stubReceipt('rgi_collect', org)
}

process CPSS_PREFILTER {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('13_cpss.py')} prefilter ${pair(meta)}
    ${receipt('cpss_prefilter', meta.id)}
    """

    stub:
    stubReceipt('cpss_prefilter', meta.id)
}

process CPSS_CHUNK {
    tag "${meta.id} ${chunk}"
    label 'cpss'

    input:
    tuple val(meta), val(chunk), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('13_cpss.py')} run ${pair(meta)} --chunk ${chunk} --threads ${task.cpus}
    ${receipt('cpss_chunk', "${meta.id} ${chunk}")}
    """

    stub:
    stubReceipt('cpss_chunk', "${meta.id} ${chunk}")
}

process CPSS_SELECT {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('13_cpss.py')} select ${pair(meta)}
    ${receipt('candidates', meta.id)}
    """

    stub:
    stubReceipt('candidates', meta.id)
}

process PREVALENCE {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('10_prevalence.py')} ${pair(meta)}
    ${receipt('prevalence', meta.id)}
    """

    stub:
    stubReceipt('prevalence', meta.id)
}

process MDA {
    tag "${meta.id}"
    label 'mda'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('12_mda.py')} ${pair(meta)}
    ${receipt('mda', meta.id)}
    """

    stub:
    stubReceipt('mda', meta.id)
}

process CARD_LAYER {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('09_card_layer.py')} ${pair(meta)}
    ${receipt('card_layer', meta.id)}
    """

    stub:
    stubReceipt('card_layer', meta.id)
}

process PYSEER_PREP {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('run_pyseer.sh'), emit: script

    script:
    """
    ${py('14_pyseer.py')} prep ${pair(meta)} --threads ${task.cpus}
    cp "\$(${py('config_path.py')} pyseer_dir ${pair(meta)})/run_pyseer.sh" run_pyseer.sh
    """

    stub:
    """
    echo 'exit 0' > run_pyseer.sh
    """
}

process PYSEER_LMM {
    tag "${meta.id}"
    label 'tools'

    input:
    tuple val(meta), path(script)

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    bash ${script}
    ${receipt('pyseer_lmm', meta.id)}
    """

    stub:
    stubReceipt('pyseer_lmm', meta.id)
}

process PYSEER_POST {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('14_pyseer.py')} post ${pair(meta)}
    ${receipt('pyseer', meta.id)}
    """

    stub:
    stubReceipt('pyseer', meta.id)
}

process GRADING {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('14b_grading.py')} ${pair(meta)}
    ${receipt('grading', meta.id)}
    """

    stub:
    stubReceipt('grading', meta.id)
}

process LP_CHUNK {
    tag "${meta.id} f${fold} c${chunk}"
    label 'permutation'

    input:
    tuple val(meta), val(fold), val(chunk), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('12b_label_permutation.py')} run ${pair(meta)} --fold ${fold} --chunk ${chunk} \\
        --threads ${task.cpus}
    ${receipt('label_permutation_chunk', "${meta.id} ${fold} ${chunk}")}
    """

    stub:
    stubReceipt('label_permutation_chunk', "${meta.id} ${fold} ${chunk}")
}

process LP_METRICS {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(meta), path('receipt.json'), emit: done

    script:
    """
    ${py('12b_label_permutation.py')} metrics ${pair(meta)}
    ${receipt('label_permutation', meta.id)}
    """

    stub:
    stubReceipt('label_permutation', meta.id)
}

process LP_ACROSS {
    label 'light'

    input:
    path deps, stageAs: 'dep*.json'

    output:
    path 'receipt.json', emit: done

    script:
    """
    ${py('12b_label_permutation.py')} across
    ${receipt('label_permutation_across', 'panel')}
    """

    stub:
    stubReceipt('label_permutation_across', 'panel')
}

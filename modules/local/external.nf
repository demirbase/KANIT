// Comparison with genotype-based prediction, per organism (16).

include { py; receipt; stubReceipt } from './common'

process EXTERNAL_PREP {
    tag "${org}"
    label 'light'

    input:
    tuple val(org), val(antibiotics), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('run_external.sh'), emit: script
    path 'receipt.json', emit: receipt

    script:
    """
    ${py('16_external.py')} prep --organism ${org} --threads ${task.cpus}
    cp "\$(${py('config_path.py')} external_dir --organism ${org})/run_external.sh" run_external.sh
    ${receipt(task, 'external_prep', org)}
    """

    stub:
    """
    echo 'exit 0' > run_external.sh
    ${stubReceipt(task, 'external_prep', org)}
    """
}

process EXTERNAL_RUN {
    tag "${org} ${shard}/${n}"
    label 'tools'

    input:
    tuple val(org), path(script), val(shard), val(n)

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    bash ${script} ${shard} ${n}
    ${receipt(task, 'external_tools', "${org} ${shard}/${n}")}
    """

    stub:
    stubReceipt(task, 'external_tools', "${org} ${shard}/${n}")
}

process EXTERNAL_COLLECT {
    tag "${org}"
    label 'light'

    input:
    tuple val(org), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('16_external.py')} collect --organism ${org}
    ${receipt(task, 'external_collect', org)}
    """

    stub:
    stubReceipt(task, 'external_collect', org)
}

process EXTERNAL_COMPARE {
    tag "${org}"
    label 'light'

    input:
    tuple val(org), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('16_external.py')} compare --organism ${org}
    ${receipt(task, 'external_compare', org)}
    """

    stub:
    stubReceipt(task, 'external_compare', org)
}

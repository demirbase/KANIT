// Comparison with genotype-based prediction, per organism (16).

include { py; receipt; stubReceipt } from './common'

process EXTERNAL_PREP {
    tag "${org}"
    label 'light'

    input:
    tuple val(org), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('run_external.sh'), emit: script

    script:
    """
    ${py('16_external.py')} prep --organism ${org} --threads ${task.cpus}
    cp "\$(${py('config_path.py')} external_dir --organism ${org})/run_external.sh" run_external.sh
    """

    stub:
    """
    echo 'exit 0' > run_external.sh
    """
}

process EXTERNAL_RUN {
    tag "${org}"
    label 'tools'

    input:
    tuple val(org), path(script)

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    bash ${script}
    ${receipt('external_tools', org)}
    """

    stub:
    stubReceipt('external_tools', org)
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
    ${receipt('external_collect', org)}
    """

    stub:
    stubReceipt('external_collect', org)
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
    ${receipt('external_compare', org)}
    """

    stub:
    stubReceipt('external_compare', org)
}

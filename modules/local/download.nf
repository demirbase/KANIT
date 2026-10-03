// Genomes and phenotypes from BV-BRC (00a, 00). Needs internet: the DOWNLOAD entry
// runs it on a login node (F3.3 settles the snapshot rules).

include { py; receipt; stubReceipt } from './common'

process DOWNLOAD_BVBRC {
    tag "${org}"
    label 'internet'

    input:
    val org

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('00a_download_bvbrc.py')} --organism ${org}
    ${receipt('download', org)}
    """

    stub:
    stubReceipt('download', org)
}

process PREPARE_METADATA {
    tag "${org}"
    label 'light'

    input:
    tuple val(org), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('00_prepare_metadata.py')} --organism ${org}
    ${receipt('metadata', org)}
    """

    stub:
    stubReceipt('metadata', org)
}

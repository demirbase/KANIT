// The data snapshot of an organism from BV-BRC (00a): genomes, phenotypes and
// assemblies, frozen with their checksums. Needs internet: the DOWNLOAD entry runs it
// on a login node. With a frozen snapshot it only verifies the checksums.

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
    ${py('00a_download_bvbrc.py')} all --organism ${org}
    ${receipt(task, 'download', org)}
    """

    stub:
    stubReceipt(task, 'download', org)
}

// NCBI context of an organism's candidate unitigs (18). Not evidence.

include { py; receipt; stubReceipt } from './common'

// needs internet: the CONTEXT entry runs it on a login node
process CONTEXT_QUERY {
    tag "${org}"
    label 'internet'
    maxForks 1                      // NCBI: one search at a time

    input:
    val org

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('18_ncbi_context.py')} query --organism ${org}
    ${receipt('context_query', org)}
    """

    stub:
    stubReceipt('context_query', org)
}

process CONTEXT_BUILD {
    tag "${org}"
    label 'light'

    input:
    val org

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    def allow = params.allow_missing_context ? '--allow-incomplete' : ''
    """
    ${py('18_ncbi_context.py')} build --organism ${org} ${allow}
    ${receipt('context_build', org)}
    """

    stub:
    stubReceipt('context_build', org)
}

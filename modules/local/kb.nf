// The knowledge base (build_kb.py), the hypothesis tests (hypotheses.py) and the release
// package of the knowledge base (kb_package.py).

include { py; receipt; stubReceipt } from './common'

process BUILD_KB {
    label 'kb'

    input:
    path deps, stageAs: 'dep*.json'

    output:
    path 'receipt.json', emit: done

    script:
    def allow = params.allow_missing_context ? '--allow-missing-context' : ''
    // the licence and the DOI of a release; a build without them packages as a draft
    def release = (params.kb_license ? " --license '${params.kb_license}'" : '') +
                  (params.kb_doi ? " --doi '${params.kb_doi}'" : '')
    """
    ${py('build_kb.py')} ${allow}${release}
    ${receipt(task, 'knowledge_base', 'kanit')}
    """

    stub:
    stubReceipt(task, 'knowledge_base', 'kanit')
}

process HYPOTHESES {
    label 'light'

    input:
    path deps, stageAs: 'dep*.json'

    output:
    path 'receipt.json', emit: done

    script:
    """
    ${py('hypotheses.py')}
    ${receipt(task, 'hypotheses', 'kanit')}
    """

    stub:
    stubReceipt(task, 'hypotheses', 'kanit')
}

process KB_PACKAGE {
    label 'light'

    input:
    path deps, stageAs: 'dep*.json'

    output:
    path 'receipt.json', emit: done

    script:
    """
    ${py('kb_package.py')}
    ${receipt(task, 'kb_package', 'kanit')}
    """

    stub:
    stubReceipt(task, 'kb_package', 'kanit')
}

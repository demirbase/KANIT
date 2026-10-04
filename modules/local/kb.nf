// The knowledge base (build_kb.py) and the hypothesis tests (hypotheses.py).

include { py; receipt; stubReceipt } from './common'

process BUILD_KB {
    label 'kb'

    input:
    path deps, stageAs: 'dep*.json'

    output:
    path 'receipt.json', emit: done

    script:
    def allow = params.allow_missing_context ? '--allow-missing-context' : ''
    """
    ${py('build_kb.py')} ${allow}
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

// Genomes of an organism: quality control (02d with CheckM2 and QUAST), lineages
// (02c), the panel of every organism (02e) and the unitig store (03u).

include { py; receipt; stubReceipt } from './common'

process QC_PREP {
    tag "${org}"
    label 'light'

    input:
    val org

    output:
    tuple val(org), path('qc_paths.sh'), emit: paths
    path 'receipt.json', emit: receipt

    script:
    """
    ${py('02d_genome_qc.py')} --mode prep --organism ${org}
    cp "\$(${py('config_path.py')} genome_qc_dir --organism ${org})/02d_qc_paths_${org}.sh" qc_paths.sh
    ${receipt(task, 'qc_prep', org)}
    """

    stub:
    """
    echo 'GENOMES_DIR=stub' > qc_paths.sh
    ${stubReceipt(task, 'qc_prep', org)}
    """
}

process CHECKM2 {
    tag "${org}"
    label 'checkm2'

    input:
    tuple val(org), path(paths)

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    source ${paths}
    checkm2 predict --input "\$GENOMES_DIR" -x fna --output-directory "\$CHECKM2_OUT" \\
        --threads ${task.cpus} --force
    checkm2 --version > "\$CHECKM2_OUT/version.txt"
    ${receipt(task, 'checkm2', org, [checkm2: 'checkm2 --version'])}
    """

    stub:
    stubReceipt(task, 'checkm2', org)
}

process QUAST {
    tag "${org}"
    label 'tools'

    input:
    tuple val(org), path(paths)

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    source ${paths}
    quast.py "\$GENOMES_DIR"/*.fna -o "\$QUAST_OUT" --threads ${task.cpus} --no-plots
    quast.py --version > "\$QUAST_OUT/version.txt"
    ${receipt(task, 'quast', org, [quast: 'quast.py --version'])}
    """

    stub:
    stubReceipt(task, 'quast', org)
}

process QC_POST {
    tag "${org}"
    label 'light'

    input:
    tuple val(org), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('02d_genome_qc.py')} --mode post --organism ${org}
    ${receipt(task, 'genome_qc', org)}
    """

    stub:
    stubReceipt(task, 'genome_qc', org)
}

process LINEAGE {
    tag "${org}"
    label 'lineage'

    input:
    val org

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('02c_lineage_poppunk.py')} --organism ${org} --threads ${task.cpus}
    ${receipt(task, 'lineage', org)}
    """

    stub:
    stubReceipt(task, 'lineage', org)
}

process PANEL {
    label 'light'

    input:
    val orgs
    path deps, stageAs: 'dep*.json'

    output:
    path 'panel_decisions.csv', emit: decisions
    path 'receipt.json', emit: receipt

    script:
    """
    ${py('02e_panel.py')} --organisms ${orgs.join(' ')}
    cp "\$(${py('config_path.py')} panel_dir)/panel_decisions.csv" panel_decisions.csv
    ${receipt(task, 'panel', orgs.join(','))}
    """

    stub:
    """
    cp ${params.stub_panel} panel_decisions.csv
    ${stubReceipt(task, 'panel', orgs.join(','))}
    """
}

process STORE {
    tag "${org}"
    label 'store'

    input:
    tuple val(org), val(antibiotics), path(deps, stageAs: 'dep*.json')

    output:
    tuple val(org), path('receipt.json'), emit: done

    script:
    """
    ${py('03u_unitig_matrix.py')} --organism ${org} --build-store --threads ${task.cpus}
    ${receipt(task, 'unitig_store', org)}
    """

    stub:
    stubReceipt(task, 'unitig_store', org)
}

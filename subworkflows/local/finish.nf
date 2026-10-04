// The end of every entry workflow: the completeness report, the checksums of the outputs,
// the verified backup of what changed (unless --backup_remote none or unset), SLURM's
// accounting of the run, the resource report and, last, the completeness gate, which
// fails an incomplete run. Everything is written in the run directory (params.trace_dir).

include { COMPLETENESS; RUN_OUTPUTS; BACKUP_PACK; BACKUP_UPLOAD; SACCT_DUMP; RUN_RESOURCES;
          COMPLETENESS_GATE } from '../../modules/local/run'

workflow FINISH {
    take:
    receipts                        // the receipts of the entry workflow's last tasks, collected
    entry                           // main, DOWNLOAD, CONTEXT or KB

    main:
    COMPLETENESS(receipts, entry)
    RUN_OUTPUTS(COMPLETENESS.out.done)
    def last = RUN_OUTPUTS.out.done
    if (params.backup_remote && params.backup_remote != 'none') {
        BACKUP_PACK(last)
        BACKUP_UPLOAD(BACKUP_PACK.out.done)
        last = BACKUP_UPLOAD.out.done
    }
    SACCT_DUMP(last)
    RUN_RESOURCES(SACCT_DUMP.out.done)
    COMPLETENESS_GATE(RUN_RESOURCES.out.done)

    emit:
    done = COMPLETENESS_GATE.out.done
}

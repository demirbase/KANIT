// The end of every entry workflow: the checksums of the outputs, the verified backup of
// what changed (unless --backup_remote none or unset), SLURM's accounting of the run
// and the resource report, all in the run directory (params.trace_dir).

include { RUN_OUTPUTS; BACKUP_PACK; BACKUP_UPLOAD; SACCT_DUMP; RUN_RESOURCES } from '../../modules/local/run'

workflow FINISH {
    take:
    receipts                        // the receipts of the entry workflow's last tasks, collected

    main:
    RUN_OUTPUTS(receipts)
    def last = RUN_OUTPUTS.out.done
    if (params.backup_remote && params.backup_remote != 'none') {
        BACKUP_PACK(last)
        BACKUP_UPLOAD(BACKUP_PACK.out.done)
        last = BACKUP_UPLOAD.out.done
    }
    SACCT_DUMP(last)
    RUN_RESOURCES(SACCT_DUMP.out.done)

    emit:
    done = RUN_RESOURCES.out.done
}

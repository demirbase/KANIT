#!/bin/bash
# The Nextflow head job of a run of many days on TRUBA (conf/truba.config). The login node
# allows a process 300 s of CPU (360 s hard), which a head spends in a day of many tasks,
# and a resumed head spends more the more tasks it finds cached; every partition ends a job
# after 3 days. This job runs the head on a barbun node. Two hours before its limit it stops
# the head with SIGKILL, so that the head's jobs keep running, cancels the jobs the head left
# waiting and submits itself again, to start once the running ones have ended and resume.
# A head on a compute node has no internet: run with --backup_remote none and back up from
# the login node afterwards.
#
#   cd $AMR_WORK && sbatch [--dependency=afterany:<job ids>] $AMR_HOME/conf/truba_head.sbatch \
#       <log file> <arguments of nextflow run main.nf, with -resume <session>>
#
#SBATCH --job-name=kanit_head
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=3-00:00:00
#SBATCH --exclude=barbun45,barbun46,barbun89,barbun118
#SBATCH --output=kanit_head_%j.log

set -u
OUT=$1
shift
STOP_AFTER=${KANIT_HEAD_HOURS:-70}          # hours before the head is stopped and resubmitted

cd "$AMR_HOME" || exit 1
export NXF_OPTS='-Xms1g -Xmx4g -XX:+UseSerialGC -XX:-UsePerfData'
echo "$(date '+%F %T') head job $SLURM_JOB_ID on $(hostname): nextflow run main.nf $*" >> "$OUT"
nextflow run main.nf "$@" >> "$OUT" 2>&1 < /dev/null &
NF=$!

for (( t = 0; t < STOP_AFTER * 3600; t += 60 )); do
    if ! kill -0 "$NF" 2>/dev/null; then     # the run ended, completed or failed
        wait "$NF"
        exit $?
    fi
    sleep 60
done

kill -9 "$NF"                                 # its jobs keep running
sleep 10
ours() { squeue -u "$USER" -h -t "$1" -o '%i %j' | awk '$2 ~ /^kanit_/ && $2 != "kanit_head" {print $1}'; }
pending=$(ours PD)
if [ -n "$pending" ]; then scancel $pending; fi
running=$(ours R,CG | paste -sd: -)
cd "$AMR_WORK" || exit 1
sbatch ${running:+--dependency=afterany:$running} "$AMR_HOME/conf/truba_head.sbatch" "$OUT" "$@"
echo "$(date '+%F %T') head stopped after ${STOP_AFTER} h; resubmitted after: ${running:-nothing}" >> "$OUT"

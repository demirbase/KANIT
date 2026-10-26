#!/bin/bash
# The Nextflow head job of a run of many days on TRUBA (conf/truba.config). The login node
# allows a process 300 s of CPU (360 s hard), which a head spends in a day of many tasks,
# and a resumed head spends more the more tasks it finds cached. This job runs the head on a barbun
# node, or in the partition given with -p (then also give -t and KANIT_HEAD_HOURS). It asks for two
# days, not the partition's three: a job that asks for the whole limit cannot be backfilled into the
# gaps the scheduler keeps for higher-priority jobs, and it then waits for hours with nodes idle
# (2026-10-15). A head is a babysitter and uses one core of the many it must ask for, but no
# partition takes a smaller job. Two days also outlast the longest task, the 36 hours of a large
# organism's units: a head that ends while its own jobs run leaves the next one waiting on them,
# and nothing is submitted meanwhile, so such a handover empties the queue (2026-10-18). It
# resubmits itself to the partition it is running in. Two hours before its limit it stops the head
# with SIGKILL, so that the head's jobs keep running, cancels the jobs the head left waiting and
# submits itself again, to start once the running ones have ended and resume.
# A head on a compute node has no internet: run with --backup_remote none and back up from
# the login node afterwards. While it waits for a node, the head may run on the login node
# (conf/truba_login_head.sh), which stops when this job runs; this job then cancels the jobs
# left waiting and starts once the running ones have ended. A second head job that finds one
# running ends at once.
#
#   cd $AMR_WORK && sbatch [--dependency=afterany:<job ids>] $AMR_HOME/conf/truba_head.sbatch \
#       <log file> <arguments of nextflow run main.nf, with -resume <session>>
#
#SBATCH --job-name=kanit_head
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2-00:00:00
#SBATCH --exclude=barbun45,barbun46,barbun89,barbun106,barbun118,hamsi20,hamsi29
#SBATCH --output=kanit_head_%j.log

set -u
OUT=$1
shift
STOP_AFTER=${KANIT_HEAD_HOURS:-22}          # hours before the head is stopped and resubmitted
POLL=${KANIT_POLL:-60}                      # seconds between checks
ours() { squeue -u "$USER" -h -t "$1" -o '%i %j' | awk '$2 ~ /^kanit_/ && $2 != "kanit_head" {print $1}'; }

# One head at a time. Several head jobs may be queued at once, in different partitions, so that
# the first to find a node wins and the others end: which partition frees up cannot be predicted,
# and a job queued in one of them is free. Asking the queue whether another head runs is not
# enough, though: two head jobs started in the same scheduling cycle can each find the other
# still pending, and two Nextflow heads on one work directory submit every task twice and write
# over each other's cache. mkdir is atomic on the shared file system, so exactly one of them
# wins. The queue is asked only about a lock that is already held, to tell a running owner from
# one that is gone: a head killed at its time limit leaves the directory behind.
# Taking it over needs proof that its owner has ended, from sacct: the owner's id is written
# just after the directory is made, and a job that reads the lock in that gap finds no owner at
# all, while a job that has just started is not yet in squeue's running jobs. Both look stale and
# are not. An unclear answer therefore ends this head, which costs a resubmission, where a wrong
# takeover would cost the cache. sacct may be missing on a compute node, and then no head can take
# the lock over: the log line says which directory to remove.
LOCK=$AMR_WORK/kanit_head.lock
take_lock() { mkdir "$LOCK" 2>/dev/null && echo "${SLURM_JOB_ID:-$$}" > "$LOCK/owner"; }
if ! take_lock; then
    owner=$(cat "$LOCK/owner" 2>/dev/null)
    state=$(sacct -n -X -j "${owner:-0}" -o State 2>/dev/null | head -1 | awk '{print $1}')
    case "$state" in
    COMPLETED|FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|BOOT_FAIL|DEADLINE|REVOKED)
        echo "$(date '+%F %T') head job ${SLURM_JOB_ID:-}: the lock of $owner is stale ($state)" >> "$OUT"
        rm -rf "$LOCK"
        take_lock || { echo "$(date '+%F %T') head job ${SLURM_JOB_ID:-}: another head took the lock; this one ends" >> "$OUT"; exit 0; } ;;
    *)
        echo "$(date '+%F %T') head job ${SLURM_JOB_ID:-}: head job ${owner:-unknown} holds the lock (${state:-no state}); this one ends. If no head runs, remove $LOCK" >> "$OUT"
        exit 0 ;;
    esac
fi
trap 'rm -rf "$LOCK"' EXIT
while [ -e "$AMR_WORK/kanit_login_head" ]; do sleep "$POLL"; done   # the login head stops
pending=$(ours PD)
if [ -n "$pending" ]; then scancel $pending; fi
while [ -n "$(ours R,CG)" ]; do sleep "$POLL"; done                    # its running jobs end

cd "$AMR_HOME" || exit 1
export NXF_OPTS='-Xms1g -Xmx4g -XX:+UseSerialGC -XX:-UsePerfData'
echo "$(date '+%F %T') head job $SLURM_JOB_ID on $(hostname): nextflow run main.nf $*" >> "$OUT"
nextflow run main.nf "$@" >> "$OUT" 2>&1 < /dev/null &
NF=$!

for (( t = 0; t < STOP_AFTER * 3600; t += POLL )); do
    if ! kill -0 "$NF" 2>/dev/null; then     # the run ended, completed or failed
        wait "$NF"
        exit $?
    fi
    sleep "$POLL"
done

kill -9 "$NF"                                 # its jobs keep running
pkill -KILL -u "$USER" -f 'run main.nf'
sleep 10
pending=$(ours PD)
if [ -n "$pending" ]; then scancel $pending; fi
running=$(ours R,CG | paste -sd: -)
cd "$AMR_WORK" || exit 1
# the next head stays in the partition where this one found a place: a partition that takes the job
# within minutes is worth more than a longer one that leaves it queued for hours (2026-10-15). Each
# partition has a smallest job it accepts, and a resubmission that asks for fewer cores than that is
# refused, so the core count travels with the partition: barbun 20, hamsi 56, orfoz 56. TRUBA's
# queue table says hamsi takes 28, but the queue itself refuses such a job and the message it
# prints asks for 56, a whole node, which is what conf/truba.config has always given it
# (2026-10-20). orfoz is the largest partition, 192 nodes of 112 cores, and the only one with
# room when barbun and hamsi are full (2026-10-20). Every partition but debug gets the two days
# that outlast a large organism's units; asking for one day instead bought no place in the
# queue (squeue --start planned none) and would stall the run at the handover.
case "${SLURM_JOB_PARTITION:-barbun}" in
    debug) PART=debug;  LIMIT=4:00:00;     HOURS=3;  CPUS=20;;
    hamsi) PART=hamsi;  LIMIT=2-00:00:00;  HOURS=46; CPUS=56;;
    orfoz) PART=orfoz;  LIMIT=2-00:00:00;  HOURS=46; CPUS=56;;
    *)     PART=barbun; LIMIT=2-00:00:00;  HOURS=46; CPUS=20;;
esac
sbatch -p "$PART" -t "$LIMIT" -c "$CPUS" --export=ALL,KANIT_HEAD_HOURS="$HOURS" \
    ${running:+--dependency=afterany:$running} "$AMR_HOME/conf/truba_head.sbatch" "$OUT" "$@"
echo "$(date '+%F %T') head stopped after ${STOP_AFTER} h; resubmitted to $PART after: ${running:-nothing}" >> "$OUT"

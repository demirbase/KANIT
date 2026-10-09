#!/bin/bash
# The Nextflow head on the login node while a head job (conf/truba_head.sbatch) waits for a
# node. The login node allows a process 300 s of CPU (360 s hard), so the head is stopped as
# soon as a head job runs: SIGKILL (its jobs keep running), its waiting jobs cancelled and
# the marker the head job waits for removed; the head job resumes once the running jobs
# have ended. If the head dies first, the marker goes too.
#
#   setsid nohup bash $AMR_HOME/conf/truba_login_head.sh <log file> <arguments of nextflow
#       run main.nf, with -resume <session>> > /dev/null 2>&1 < /dev/null &

set -u
OUT=$1
shift
MARK=$AMR_WORK/kanit_login_head
POLL=${KANIT_POLL:-60}                      # seconds between checks

cd "$AMR_HOME" || exit 1
ulimit -St 360 2>/dev/null
export NXF_OPTS='-Xms1g -Xmx4g -XX:+UseSerialGC -XX:-UsePerfData'
echo "$(date '+%F %T') head on the login node: nextflow run main.nf $*" >> "$OUT"
nextflow run main.nf "$@" >> "$OUT" 2>&1 < /dev/null &
NF=$!
echo "$NF" > "$MARK"

while kill -0 "$NF" 2>/dev/null; do
    if squeue -u "$USER" -h -n kanit_head -t R -o '%i' | grep -q .; then
        kill -9 "$NF"
        pkill -KILL -u "$USER" -f 'run main.nf'
        sleep 10
        squeue -u "$USER" -h -t PD -o '%i %j' |
            awk '$2 ~ /^kanit_/ && $2 != "kanit_head" {print $1}' | xargs -r scancel
        echo "$(date '+%F %T') login head stopped: a head job runs" >> "$OUT"
        break
    fi
    sleep "$POLL"
done
rm -f "$MARK"

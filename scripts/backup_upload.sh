#!/usr/bin/env bash
# Upload and verification of a packed backup (scripts/backup.py pack; plan F3.5).
# Runs on a node with internet: on TRUBA the login node, outside the containers.
#
#   backup_upload.sh STAGE REMOTE RUN LEDGER
#
# One archive per rclone call: a short process is not killed for running long, and a
# failure costs one archive. Then three checks against REMOTE/RUN:
#   1. the MD5 of every archive and of the manifest (rclone check, scoped to them);
#   2. no name twice (Google Drive allows two files with one name, and a retried upload
#      can make them);
#   3. the number of files.
# Only when all three pass are the manifest's rows added to the ledger, the ledger
# copied to REMOTE and the staged archives deleted. Otherwise it exits non-zero and
# leaves everything in place: a backup is never reported done without being checked.
# rclone exit 7 is Google Drive's daily upload cap: run the stage again the next day.
set -uo pipefail
[[ $# -eq 4 ]] || { echo "usage: $0 STAGE REMOTE RUN LEDGER" >&2; exit 2; }
STAGE=$1; REMOTE=${2%/}; RUN=$3; LEDGER=$4
RCLONE="${RCLONE:-rclone}"
MANIFEST="$STAGE/backup_manifest.tsv"
[[ -f "$MANIFEST" ]] || { echo "ERROR: no $MANIFEST (run backup.py pack first)" >&2; exit 2; }

ARCHIVES=()
while IFS=$'\t' read -r _ _ _ _ _ archive _; do
    ARCHIVES+=("$archive")
done < <(tail -n +2 "$MANIFEST")
if [[ ${#ARCHIVES[@]} -eq 0 ]]; then
    echo "BACKUP UPLOAD — nothing changed since the last verified backup"
    exit 0
fi

DEST="$REMOTE/$RUN"
OPTS=(--transfers 1 --checkers 2 --buffer-size 0 --use-mmap --drive-chunk-size 32M
      --drive-stop-on-upload-limit --retries 5 --low-level-retries 20
      --stats 2m --stats-one-line)
failed=()
for a in "${ARCHIVES[@]}" backup_manifest.tsv; do
    "$RCLONE" copy "$STAGE/$a" "$DEST/" "${OPTS[@]}"; rc=$?
    case $rc in
        0) echo "  up $a" ;;
        7) echo "ERROR: the daily upload cap was reached at $a; run the stage again tomorrow" >&2
           exit 7 ;;
        *) echo "  FAILED $a (rclone exit $rc)" >&2; failed+=("$a") ;;
    esac
done
if [[ ${#failed[@]} -gt 0 ]]; then
    echo "ERROR: ${#failed[@]} upload(s) failed: ${failed[*]}" >&2
    exit 1
fi

LIST=$(mktemp "$STAGE/.include.XXXXXX")
printf '%s\n' "${ARCHIVES[@]}" backup_manifest.tsv > "$LIST"
if ! "$RCLONE" check "$STAGE" "$DEST" --one-way --include-from "$LIST" --checkers 2; then
    rm -f "$LIST"
    echo "ERROR: the files on $DEST do not match the staged ones" >&2
    exit 1
fi
rm -f "$LIST"

NAMES=$("$RCLONE" lsf "$DEST" --files-only) || { echo "ERROR: cannot list $DEST" >&2; exit 1; }
DUP=$(printf '%s\n' "$NAMES" | sort | uniq -d)
if [[ -n "$DUP" ]]; then
    echo "ERROR: names on $DEST twice: $DUP (rclone dedupe --by-hash)" >&2
    exit 1
fi
N=$(printf '%s\n' "$NAMES" | grep -c .)
EXPECTED=$(( ${#ARCHIVES[@]} + 1 ))
if [[ "$N" -ne "$EXPECTED" ]]; then
    echo "ERROR: $N files on $DEST, expected $EXPECTED" >&2
    exit 1
fi

NOW=$(date -u +%Y-%m-%dT%H:%M:%SZ)
if [[ ! -f "$LEDGER" ]]; then
    mkdir -p "$(dirname "$LEDGER")"
    printf 'unit\tmode\tn_files\tbytes\tfingerprint\tarchive\tarchive_bytes\tmd5\tremote\trun\tverified_at\n' > "$LEDGER"
fi
awk -F'\t' -v OFS='\t' -v dest="$DEST" -v run="$RUN" -v now="$NOW" \
    'NR > 1 {print $0, dest, run, now}' "$MANIFEST" >> "$LEDGER"
"$RCLONE" copy "$LEDGER" "$REMOTE/" "${OPTS[@]}" \
    || { echo "ERROR: the ledger could not be copied to $REMOTE" >&2; exit 1; }
for a in "${ARCHIVES[@]}"; do
    rm -f "$STAGE/$a"
done
echo "BACKUP UPLOAD — ${#ARCHIVES[@]} archive(s) verified on $DEST"

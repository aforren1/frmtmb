#!/usr/bin/env bash
# Replays .github/scripts/skip-gate.sh on the testthat output of real
# CI runs, downloaded with `gh run view <id> --log` into
# dev/smallfix-ci-log/<name>/run.log. For each job, the "Skipped tests"
# block that the job printed, and testthat's summary line, are written
# as a testthat.Rout under a directory laid out as the job's check or
# covr directory, and the gate runs on it with the allow-list given.
#   bash dev/smallfix-gate-replay.sh <run.log> <out dir> <job> [allowed..]
# Prints the gate's output and exit status.
set -u
log=$1; out=$2; job=$3; shift 3
d="$out/$(echo "$job" | tr -c 'A-Za-z0-9.\n' '_')"
rm -rf "$d"; mkdir -p "$d/frmtmb.Rcheck/tests"
# a job's lines are "<job>\t<step>\t<timestamp> <text>"; keep the text
# of every line from the skip header to the next testthat header, and
# the last summary line
awk -F'\t' -v j="$job" '$1 == j { sub(/^[0-9TZ:.-]+ /, "", $3); print $3 }' \
  "$log" > "$d/job.txt"
awk '/══ Skipped tests/ { on = 1; print; next }
     on && /^══ / { on = 0 }
     on { print }' "$d/job.txt" > "$d/frmtmb.Rcheck/tests/testthat.Rout"
grep -E '^\[ FAIL [0-9]+ \| WARN' "$d/job.txt" | tail -1 \
  >> "$d/frmtmb.Rcheck/tests/testthat.Rout"
rm -f "$d/job.txt"
echo "-- $job: $(grep -c '' "$d/frmtmb.Rcheck/tests/testthat.Rout") lines" \
  "replayed; allowed: ${*:-none}"
bash .github/scripts/skip-gate.sh "$d" "$@"
echo "-- exit $?"

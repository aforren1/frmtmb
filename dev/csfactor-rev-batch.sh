#!/bin/sh
# Reviewer batch: one test file per R process, three at a time, every
# RESULT line appended to one log. A file that produces no RESULT line is
# reported as NO-RESULT so an aborted process cannot read as a clean pass.
#   sh dev/csfactor-rev-batch.sh <lib> <log> <file> [<file> ...]
LIB="$1"; LOG="$2"; shift 2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
WT=C:/Users/adf44/source/r/frmtmb-wt-csfactor
: > "$LOG"
n=0
for f in "$@"; do
  (
    out=$(NOT_CRAN=true "$R" "$WT/dev/csfactor-rev-run.R" "$LIB" "$f" 2>&1)
    if echo "$out" | grep -q '^RESULT'; then
      echo "$out" | grep -E '^(RESULT|  BAD|  SKIP|     )' >> "$LOG"
    else
      echo "RESULT $f NO-RESULT-LINE" >> "$LOG"
      echo "$out" | tail -25 >> "$LOG"
    fi
  ) &
  n=$((n + 1))
  if [ $((n % 3)) -eq 0 ]; then wait; fi
done
wait
echo "BATCH files requested: $#" >> "$LOG"
echo "BATCH RESULT lines: $(grep -c '^RESULT' "$LOG")" >> "$LOG"

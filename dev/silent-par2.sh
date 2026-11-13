#!/usr/bin/env bash
# Runs "arm pkg testpath" job lines from a file, at most 3 R processes.
# Usage: bash dev/silent-par.sh <jobs file> <log dir>
jobs_file=$1
logdir=$2
mkdir -p "$logdir"
RS="${RS:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
cd /c/Users/adf44/source/r/frmtmb-wt-silent
while read -r arm pkg path; do
  [ -z "$arm" ] && continue
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 2; done
  name=$(basename "$path" .R)
  "$RS" dev/silent-run-test.R "$arm" "$pkg" "$path" > "$logdir/$arm-$pkg-$name.log" 2>&1 &
done < "$jobs_file"
wait
grep -h "^RESULT\|^lib:" "$logdir"/*.log

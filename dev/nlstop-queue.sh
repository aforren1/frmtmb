#!/usr/bin/env bash
# Run a job list with at most N R processes at once (lane nlstop).
# Usage: bash dev/nlstop-queue.sh <jobs file> <N>
# Each job line: <log path>|<env assignments, space separated>|<script> <args>
# Prints QUEUE DONE when every job has ended.
jobs="$1"
n="${2:-2}"
cd /c/Users/adf44/source/r/frmtmb-wt-nlstop || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
running=0
while IFS='|' read -r log envs cmd; do
  [ -z "$log" ] && continue
  case "$log" in \#*) continue ;; esac
  mkdir -p "$(dirname "$log")"
  ( env $envs "$RS" $cmd > "$log" 2>&1; echo "JOB END $log" ) &
  running=$((running + 1))
  if [ "$running" -ge "$n" ]; then
    wait -n
    running=$((running - 1))
  fi
done < "$jobs"
wait
echo "QUEUE DONE $jobs"

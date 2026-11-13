#!/usr/bin/env bash
# Review of lane names: run queued jobs, N at a time.
# Usage: bash dev/reviews/names-review-queue.sh <N> <jobfile>
# Each job line: <log> <script.R> [args...]
cd /c/Users/adf44/source/r/frmtmb-wt-names || exit 1
N=$1; jobs=$2
n=0
while read -r log script args; do
  [ -z "$log" ] && continue
  bash dev/reviews/names-review-run.sh "$script" "$log" $args &
  n=$((n + 1))
  if [ $n -ge "$N" ]; then wait -n; n=$((n - 1)); fi
done < "$jobs"
wait
echo "QUEUE DONE $jobs"

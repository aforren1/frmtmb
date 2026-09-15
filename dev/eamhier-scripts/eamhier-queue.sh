#!/bin/sh
# Lane eamhier: run a job file, N replicates at a time.
#
# Each line of the job file is one replicate and carries every argument
# the replicate script takes, whitespace separated. For the default
# script, eamhier-rep.R, that is
#
#   <arm> <bound> <seed> <ns> <nt> <outfile>
#
# One process per replicate, so a worker never holds two tapes, and the
# throttle is an argument because this machine runs other lanes at the
# same time and a 12,000-row hierarchical DDM peaks near 2.5 GB.
#
# Run:
#   sh dev/eamhier-scripts/eamhier-queue.sh <jobfile> <parallel> <log> \
#      [<script> <args per job>]
set -e
JOBS="$1"
PAR="${2:-3}"
LOG="${3:-eamhier-queue.log}"
SCRIPT="${4:-dev/eamhier-scripts/eamhier-rep.R}"
NARG="${5:-6}"
RSCRIPT="/c/Program Files/R/R-4.6.1/bin/x64/Rscript.exe"
xargs -P "$PAR" -n "$NARG" "$RSCRIPT" --vanilla "$SCRIPT" \
  < "$JOBS" >> "$LOG" 2>&1
echo "queue done: $JOBS" >> "$LOG"

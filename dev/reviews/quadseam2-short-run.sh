#!/usr/bin/env bash
# Short check of lane quadseam2's fix round: run one R script pinned to
# one P-core mask (rtmb-pitfalls.md item 21) with the reference BLAS.
#   bash quadseam2-short-run.sh <script> <log name> <mask> [args...]
# <script> is relative to the worktree. Log:
# dev/reviews/quadseam2-short-log/<log name>.txt
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-quadseam2
scr=$1; log=$2; mask=$3; shift 3
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
export R_PARALLEL_PORT=$((11000 + RANDOM % 5000))
cd "$WT"
mkdir -p dev/reviews/quadseam2-short-log
cmd //c "dev\\reviews\\quadseam2-review-aff.cmd $mask $scr $*" \
  > "dev/reviews/quadseam2-short-log/$log.txt" 2>&1
echo "exit $? $log"

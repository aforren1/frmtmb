#!/usr/bin/env bash
# Review of lane arspat2: run one review script. Usage:
#   bash dev/reviews/arspat2-review-run.sh <script.R> <log> [base] [port]
# "base" sets REVLIB="" (rellib-r21 alone).
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-arspat2
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true
if [ "${3:-}" = "base" ]; then export REVLIB=""; fi
export R_PARALLEL_PORT="${4:-11731}"
mkdir -p "$WT/dev/reviews/arspat2-review-out"
cd "$WT" && "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$1" ${REV_ARGS:-} > "$2" 2>&1
echo "exit $?" >> "$2"

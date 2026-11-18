#!/usr/bin/env bash
# Review of lane alapse: run one R script, output to dev/reviews/alapse-review-out/<log>.
#   bash dev/reviews/alapse-review-run.sh <log> <script> [args...]
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-alapse
cd "$ROOT"
log="$1"; shift
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
mkdir -p dev/reviews/alapse-review-out
"$RS" "$@" > "dev/reviews/alapse-review-out/$log" 2>&1
echo "exit $? -> dev/reviews/alapse-review-out/$log"

#!/usr/bin/env bash
# Lane alapse: run one R script with the lane's environment, output to a
# log under dev/alapse-out/.
#   bash dev/alapse-run.sh <log name> <script> [args...]
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
: "${TMP:?TMP is unset}"
mkdir -p dev/alapse-out
"$RS" "$@" > "dev/alapse-out/$log" 2>&1
echo "exit $? -> dev/alapse-out/$log"

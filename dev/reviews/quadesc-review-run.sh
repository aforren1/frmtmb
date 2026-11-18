#!/usr/bin/env bash
# Review of lane quadesc: run one R script under the reference BLAS or an
# OpenBLAS build of dev/ciharden-openblas.sh, log to dev/reviews/<log>.txt.
#   bash dev/reviews/quadesc-review-run.sh <ref|0.3.32|0.3.26> <log> <script> [args]
set -u
VER="$1"; LOG="$2"; shift 2
WT=/c/Users/adf44/source/r/frmtmb-wt-quadesc
cd "$WT"
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
if [ "$VER" = ref ]; then R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
else R="dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"; fi
"$R" "$@" > "dev/reviews/$LOG.txt" 2>&1
echo "exit $?" >> "dev/reviews/$LOG.txt"

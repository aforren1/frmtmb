#!/usr/bin/env bash
# Lane slopetrig: run one R script with the OpenBLAS copy of R that
# dev/ciharden-openblas.sh <ver> lapack built, 4 OpenBLAS threads as
# on the ubuntu runners.
# Usage: bash dev/slopetrig-ob.sh <ver> <log> <script> [args...]
ver="$1"; log="$2"; shift 2
WT=/c/Users/adf44/source/r/frmtmb-wt-slopetrig
RS="$WT/dev/ciharden-out/Rob$ver-lapack/bin/Rscript.exe"
[ -x "$RS" ] || { echo "no $RS"; exit 1; }
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-slopetrig/dev/stan-cache
cd "$WT" || exit 1
"$RS" "$@" > "$log" 2>&1
echo "EXIT $?" >> "$log"

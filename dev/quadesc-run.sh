#!/usr/bin/env bash
# Lane quadesc: run one R script under a given R build and core type.
#   bash dev/quadesc-run.sh <ref|0.3.26|0.3.32> <coretype|default> <log> <script> [args...]
# The log lands in dev/quadesc-log/<log>.txt; the exit code is appended.
set -u
VER="$1"; CT="$2"; LOG="$3"; shift 3
WT=/c/Users/adf44/source/r/frmtmb-wt-quadesc
cd "$WT"
PATH=$(printf '%s' "$PATH" | tr ':' '\n' | grep -vi rtools | paste -sd:)
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
mkdir -p dev/quadesc-log
if [ "$VER" = ref ]; then
  R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
else
  R="dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"
fi
if [ "$CT" = default ]; then unset OPENBLAS_CORETYPE; else
  export OPENBLAS_CORETYPE="$CT"; fi
"$R" "$@" > "dev/quadesc-log/$LOG.txt" 2>&1
echo "exit $?" >> "dev/quadesc-log/$LOG.txt"

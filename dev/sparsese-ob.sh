#!/usr/bin/env bash
# Lane sparsese: test files under an OpenBLAS copy of R built by
# dev/ciharden-openblas.sh, LAPACK routed too, one process per file, two
# at a time, each file under two OpenBLAS kernels.
#   bash dev/sparsese-ob.sh <version> <file> [file ...]
# Logs: dev/sparsese-out/ob-<version>-<core>-<file>.txt
set -u
VER="$1"; shift
ROOT="C:/Users/adf44/source/r/frmtmb-wt-sparsese"
RS="$ROOT/dev/ciharden-out/Rob$VER-lapack/bin/x64/Rscript.exe"
[ -f "$RS" ] || { echo "no OpenBLAS R at $RS"; exit 1; }
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-sparsese-lib
export NOT_CRAN=true
export OPENBLAS_NUM_THREADS=4
cd "$ROOT"
for f in "$@"; do
  for core in Haswell Sandybridge; do
    while [ "$(jobs -rp | wc -l)" -ge 2 ]; do wait -n; done
    o="dev/sparsese-out/ob-$VER-$core-$(basename "$f" .R).txt"
    ( OPENBLAS_CORETYPE=$core OPENBLAS_VERBOSE=2 "$RS" \
        dev/release/run-tests.R frmtmb "$f" > "$o" 2>&1 ) &
  done
done
wait
for f in "$@"; do
  for core in Haswell Sandybridge; do
    o="dev/sparsese-out/ob-$VER-$core-$(basename "$f" .R).txt"
    echo "$VER $core $(grep -a -m1 '^Core:' "$o") $(grep -a '^RESULT' "$o")"
  done
done

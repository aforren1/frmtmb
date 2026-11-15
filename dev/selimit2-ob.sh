#!/usr/bin/env bash
# Lane selimit2: test files under R with OpenBLAS (BLAS and LAPACK),
# built by dev/ciharden-openblas.sh, through dev/release/run-tests.R,
# at most N at once.
#   bash dev/selimit2-ob.sh <version> <N> <pkg:relative test file> ...
# Logs: dev/selimit2-out/ob-<version>/<pkg>--<file>.txt
set -u
VER="$1"; N="$2"; shift 2
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="$ROOT/dev/ciharden-out/Rob$VER-lapack/bin/Rscript.exe"
OUT="$ROOT/dev/selimit2-out/ob-$VER"
mkdir -p "$OUT"
export NOT_CRAN=true FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-selimit2-lib
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache" OPENBLAS_NUM_THREADS=4
for j in "$@"; do
  p=${j%%:*}; f=${j#*:}
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  "$RS" "$ROOT/dev/release/run-tests.R" "$p" "$ROOT/$f" \
    > "$OUT/$p--$(basename "$f" .R).txt" 2>&1 &
done
wait
echo "done ob $VER"

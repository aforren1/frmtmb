#!/usr/bin/env bash
# Lane code8: run test files one per process, sequentially, under an
# OpenBLAS-backed R copy built by dev/ciharden-openblas.sh.
#   bash dev/code8-ob-run.sh <0.3.26|0.3.32> <before|after> <file...>
set -u
VER="$1"; ARM="$2"; shift 2
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RS="$ROOT/dev/ciharden-out/Rob$VER-lapack/bin/x64/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
if [ "$ARM" = after ]; then
  export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-code8-lib
else
  unset FRMTMB_RELLIB
fi
mkdir -p "$ROOT/dev/code8-out/ob"
for f in "$@"; do
  o="$ROOT/dev/code8-out/ob/$VER-$ARM-$(basename "$f" .R).txt"
  "$RS" "$ROOT/dev/release/run-tests.R" frmtmb "$f" > "$o" 2>&1
  echo "$VER $ARM $(grep -a '^lib:' "$o" | cut -c1-60) | $(grep -a '^RESULT' "$o" | tail -1)"
  grep -a "DETAIL WARN\|DETAIL FAIL\|DETAIL ERR" "$o" | cut -c1-300
done

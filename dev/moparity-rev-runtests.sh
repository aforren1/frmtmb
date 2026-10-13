#!/usr/bin/env bash
# Reviewer: run test files one R process each through
# dev/release/run-tests.R with the lane's library first, logs under
# dev/moparity-rev-log/<tier>/.
#   bash dev/moparity-rev-runtests.sh <tier> <pkg:file> ...
# RSCRIPT picks the R (the OpenBLAS one for an OpenBLAS tier).
set -u
tier="$1"; shift
ROOT=C:/Users/adf44/source/r/frmtmb-wt-moparity
RS="${RSCRIPT:-/c/Program Files/R/R-4.6.1/bin/Rscript.exe}"
OUT="$ROOT/dev/moparity-rev-log/$tier"
mkdir -p "$OUT"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/moparity-rev-stan-cache"
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-moparity-lib
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
for pf in "$@"; do
  p="${pf%%:*}"; f="${pf#*:}"
  if [ "$p" = "frmtmb" ]; then d="$ROOT/tests/testthat"; else d="$ROOT/extensions/$p/tests/testthat"; fi
  ( "$RS" "$ROOT/dev/release/run-tests.R" "$p" "$d/$f" > "$OUT/$p--${f%.R}.txt" 2>&1 ) &
done
wait
for pf in "$@"; do
  p="${pf%%:*}"; f="${pf#*:}"
  o="$OUT/$p--${f%.R}.txt"
  echo "$(grep -a '^lib' "$o" | head -1) | $(grep -a '^RESULT' "$o" | tail -1)"
done

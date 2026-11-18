#!/usr/bin/env bash
# Lane alapse: test files under the OpenBLAS 0.3.26 LAPACK copy of R,
# two at a time, the lane's core first.
#   bash dev/alapse-ob.sh <file> ...
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-alapse
cd "$ROOT"
RS="$ROOT/dev/ciharden-out/Rob0.3.26-lapack/bin/x64/Rscript.exe"
RT=/c/rtools45
export PATH="$PATH:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-alapse-lib
for f in "$@"; do
  while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 2; done
  ( "$RS" dev/release/run-tests.R frmtmb "tests/testthat/$f" > "dev/alapse-out/ob026-$f.txt" 2>&1 ) &
done
wait
for f in "$@"; do grep -a "^RESULT" "dev/alapse-out/ob026-$f.txt" || echo "NO RESULT $f"; done

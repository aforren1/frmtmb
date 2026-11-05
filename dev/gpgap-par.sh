#!/usr/bin/env bash
# One test file per R process, 8 at once, through dev/gpgap-runtest.R.
#   bash dev/gpgap-par.sh <outdir> <package> [<package> ...]
# Rtools goes at the END of PATH (dev/equiv-run.sh says why).
set -u
out="$1"; shift
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export NOT_CRAN=true
mkdir -p "$out"
for p in "$@"; do
  d=extensions/$p/tests/testthat
  [ "$p" = "frmtmb" ] && d=tests/testthat
  for f in "$d"/test-*.R; do
    while [ "$(jobs -rp | wc -l)" -ge 8 ]; do sleep 1; done
    o="$out/$p--$(basename "$f" .R).txt"
    if ! grep -qa "^RESULT" "$o" 2>/dev/null; then
      ( "$RS" dev/gpgap-runtest.R "$p" "$f" > "$o" 2>&1 ) &
    fi
  done
done
wait

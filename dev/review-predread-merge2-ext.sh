#!/usr/bin/env bash
# Short check, item 7: the seven other suites on the refreshed trial merge, whose core is
# the merged one, one file per process, gated.
set -u
cd /c/Users/adf44/source/r/rev-predread-src/merge
P=/c/Users/adf44/source/r/frmtmb-wt-predread
O=$P/dev/review-predread-out/r1/merge-ext
mkdir -p "$O"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/rev-predread-src/stan-cache
LIB=C:/Users/adf44/source/r/rev-predread-lib/merge
for e in coupling eam latent learn ode spatial survival; do
  for f in extensions/frmtmb.$e/tests/testthat/test-*.R; do
    while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 2; done
    ( "$RS" "$P/dev/review-predread-runone.R" "frmtmb.$e" "$f" "$LIB" > "$O/frmtmb.$e--$(basename "$f" .R).txt" 2>&1 ) &
  done
done
wait
echo finished > "$O/../merge-ext.done"

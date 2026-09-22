#!/bin/sh
# The brms-agreement files, GATED: FRMTMB_BRMS_FIT_TESTS=true and
# NOT_CRAN=true, one R process per file, against the lane library.
# Round 2's "gated tier" was the ported suite alone and never set the
# gate for these four, which is how their breakage went unseen.
#   sh dev/shapes-gated.sh   (from the worktree root)
[ -f dev/shapes-run.R ] || { echo "run me from the worktree root"; exit 1; }
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
for f in brms-methods brms-likelihood brms-agreement brms-priors; do
  "$R" dev/shapes-run.R frmtmb "tests/testthat/test-$f.R" \
    > "dev/shapes-log/gated-$f.txt" 2>&1
  echo "exit=$? $f"
  grep -h "^RESULT" "dev/shapes-log/gated-$f.txt" || echo "NO RESULT $f"
done

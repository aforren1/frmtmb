#!/usr/bin/env bash
# Review of lane arstep: one R process at a time, in order.
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
W=C:/Users/adf44/source/r/frmtmb-wt-arstep
O=$W/dev/review-arstep-out
cd "$W" || exit 1
"$R" dev/review-arstep-brmscode.R > "$O/brmscode.txt" 2>&1
for m in control arsign nomask nosettle norestore noagree noalt noarphi nogate; do
  "$R" dev/review-arstep-mut.R "$m" > "$O/mut-$m.txt" 2>&1
done
for f in test-gp-ar-start.R test-gp-start.R test-gp-noncentered.R; do
  FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-arstep-lib NOT_CRAN=true \
    "$R" dev/arstep-run-tests.R frmtmb "tests/testthat/$f" > "$O/t-$f.txt" 2>&1
done
"$R" dev/review-arstep-freshn.R base N 0 > "$O/fresh-N-base.txt" 2>&1
echo CHAINDONE > "$O/chain-done.txt"

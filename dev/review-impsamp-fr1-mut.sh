#!/usr/bin/env bash
# Build and test mutants, at most 4 at once. Each log ends with RESULT.
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
L=review-impsamp-out/fr1/mut
mkdir -p $L
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
one() {
  m=$1
  "$RS" review-impsamp-fr1-mut.R make "$m" > $L/make-$m.log 2>&1 || return
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-impsamp-lib/$m --no-multiarch --no-docs --no-test-load C:/Users/adf44/source/r/rev-impsamp-src/fr1-$m >> $L/make-$m.log 2>&1
  "$RS" review-impsamp-fr1-mut.R test "$m" > $L/test-$m.log 2>&1
}
for m in partial_rate capped_settles abs_step no_units; do one "$m" & done
wait
for m in prior_twice fix_thres; do one "$m" & done
wait
echo done

#!/usr/bin/env bash
# One clean driver. Each mutant: source copy, install to its own library, run
# test-importance.R. The test step stops if the mutant library did not load.
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
L=review-impsamp-out/mut2
mkdir -p $L
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
one() {
  m=$1
  "$RS" review-impsamp-mutants.R make "$m" > $L/make-$m.log 2>&1 || return
  if [ "$m" = fix_h ]; then "$RS" review-impsamp-mutfix-edit.R >> $L/make-$m.log 2>&1; fi
  rm -rf "C:/Users/adf44/source/r/rev-impsamp-lib/$m/frmtmb"
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-impsamp-lib/$m --no-multiarch --no-docs --no-test-load C:/Users/adf44/source/r/rev-impsamp-src/$m >> $L/make-$m.log 2>&1
  "$RS" review-impsamp-mutants.R test "$m" test-importance.R > $L/test-$m.log 2>&1
}
n=0
for m in ${MUTS:-fix_h partial_rate no_rate no_cross frozen no_jacobian half_jacobian map_recycle no_prior_mode no_refusal one_round grad_always no_cut}; do
  one "$m" &
  n=$((n + 1))
  if [ $((n % 4)) -eq 0 ]; then wait; fi
done
wait
echo done

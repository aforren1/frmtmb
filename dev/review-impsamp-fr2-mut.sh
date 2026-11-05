#!/usr/bin/env bash
# Build, then test: no_refusal, capped_settles, partial_rate whole file; half_jacobian
# on the named off-anchor test, since its whole-file run climbs the ladder.
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
L=review-impsamp-out/fr2/mut
mkdir -p $L
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mk() {
  "$RS" review-impsamp-fr2-mut.R make "$1" > $L/make-$1.log 2>&1 &&
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=C:/Users/adf44/source/r/rev-impsamp-lib/$1 --no-multiarch --no-docs --no-test-load C:/Users/adf44/source/r/rev-impsamp-src/fr2-$1 >> $L/make-$1.log 2>&1
}
for m in partial_rate half_jacobian capped_settles no_refusal no_units theta_size no_prior_scale; do mk $m & done
wait
"$RS" review-impsamp-fr2-mut.R test partial_rate > $L/test-partial_rate.log 2>&1 &
"$RS" review-impsamp-fr2-mut.R test capped_settles > $L/test-capped_settles.log 2>&1 &
"$RS" review-impsamp-fr2-mut.R test no_refusal > $L/test-no_refusal.log 2>&1 &
"$RS" review-impsamp-fr2-mut.R test half_jacobian "the objective two standard errors off the anchor is the integral" > $L/test-half_jacobian.log 2>&1 &
wait
echo done

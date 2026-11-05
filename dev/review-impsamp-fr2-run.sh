#!/usr/bin/env bash
# All of fix round 2's checks, at most 7 processes at once. Clean libraries first:
# an install cut off at the 18:01 restart must not be read.
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp/dev || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
L=review-impsamp-out/fr2
M=$L/mut
rm -rf "$M"; mkdir -p "$M"
rm -rf /c/Users/adf44/source/r/rev-impsamp-lib /c/Users/adf44/source/r/rev-impsamp-src
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LIB=C:/Users/adf44/source/r/rev-impsamp-lib
mk() {
  "$RS" review-impsamp-fr2-mut.R make "$1" > $M/make-$1.log 2>&1 &&
  "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --library=$LIB/$1 --no-multiarch --no-docs --no-test-load C:/Users/adf44/source/r/rev-impsamp-src/fr2-$1 >> $M/make-$1.log 2>&1
}
for m in partial_rate half_jacobian capped_settles no_refusal; do mk $m & done; wait
for m in no_units theta_size no_prior_scale; do mk $m & done; wait
"$RS" review-impsamp-fr2-mut.R test partial_rate > $M/test-partial_rate.log 2>&1 &
"$RS" review-impsamp-fr2-mut.R test capped_settles > $M/test-capped_settles.log 2>&1 &
"$RS" review-impsamp-fr2-mut.R test no_refusal > $M/test-no_refusal.log 2>&1 &
"$RS" review-impsamp-fr2-mut.R test half_jacobian "the objective two standard errors off the anchor is the integral" > $M/test-half_jacobian.log 2>&1 &
REV_ARM=after "$RS" review-impsamp-fr2-natural.R same 2>&1 | grep -v "^Boundary" > $L/natural-same.txt &
REV_ARM=mutant REV_MUTLIB=$LIB/theta_size "$RS" review-impsamp-fr2-natural.R sep 2>&1 | grep -v "^Boundary" > $L/natural-sep-theta_size.txt &
REV_ARM=after "$RS" review-impsamp-fr2-natural.R search 2>&1 | grep -v "^Boundary" > $L/natural-search.txt &
wait
echo stage1
REV_ARM=after "$RS" review-impsamp-fr2-small.R units 2>&1 | grep -v "^Boundary" > $L/units-after.txt &
REV_ARM=mutant REV_MUTLIB=$LIB/no_units "$RS" review-impsamp-fr2-small.R units 2>&1 | grep -v "^Boundary" > $L/units-no_units.txt &
REV_ARM=after "$RS" review-impsamp-fr2-small.R family 2>&1 | grep -v "^Boundary" > $L/family.txt &
REV_ARM=after "$RS" review-impsamp-fr2-small.R floor 2>&1 | grep -v "^Boundary" > $L/floor.txt &
REV_ARM=after "$RS" review-impsamp-fr1-ordfix.R 2>&1 | grep -v "^Boundary" > $L/ordfix.txt &
"$RS" review-impsamp-fr2-laphigh.R 2>&1 | grep -v "^Boundary" > $L/laphigh.txt &
"$RS" review-impsamp-fr2-capped.R 2>&1 | grep -v "^Boundary" > $L/capped.txt &
wait
echo stage2
"$RS" review-impsamp-fr1-escape.R 2>&1 | grep -v "^Boundary" > $L/escape.txt &
"$RS" review-impsamp-fr1-cap.R 2>&1 | grep -v "^Boundary" > $L/cap.txt &
"$RS" review-impsamp-fr1-settled.R 2>&1 | grep -v "^Boundary" > $L/settled.txt &
"$RS" review-impsamp-fr2-cost.R 2>&1 | grep -v "^Boundary" > $L/cost.txt &
REV_ARM=after "$RS" review-impsamp-scale-se.R 2>&1 | grep -E "^==|^rounds|^SE ratio|REFUSAL|warn:" > $L/scale-se.txt &
"$RS" review-impsamp-year.R after 2>&1 | grep -v "^Boundary" > $L/year-table.txt &
wait
echo done

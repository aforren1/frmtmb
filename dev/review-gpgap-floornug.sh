#!/usr/bin/env bash
# item 7: the floor test on kernel mutants.
set -u
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
LIB="C:/Users/adf44/source/r/rev-gpgap-lib;C:/Users/adf44/source/r/rellib-r16"
for m in nug_2e6 nug_1e7 kern_102; do
  REV_MUT=$m "$RS" dev/review-gpgap-testmut.R test-gp-nugget-floor.R "$LIB" > dev/review-gpgap-out/tm-floor-$m.txt 2>&1
done

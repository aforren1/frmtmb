#!/bin/bash
# Review of lane adhess: illse2's polynomial cases on one arm, in sequence.
#   dev/adhess-rev-acc.sh <arm>
cd /c/Users/adf44/source/r/frmtmb-wt-adhess
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true REV_ARM=$1
for job in "gaussian 0" "poisson 0" "bernoulli 0" "gaussian 1"; do
  set -- $job
  "$RS" dev/adhess-rev-acc.R $1 $2 "4:8" "1:3" > dev/adhess-rev-log/acc-$REV_ARM-$1-$2.txt 2>&1
done

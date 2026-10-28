#!/bin/bash
# One test file in one R process: arm package path gated(0/1)
W=C:/Users/adf44/source/r/frmtmb-wt-reslot
OUT=$W/dev/review-reslot-out
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
arm=$1; pkg=$2; path=$3; gated=$4
log=$OUT/${5:-}t-$arm-$(basename $path .R).txt
if [ "$gated" = "1" ]; then
  export FRMTMB_BRMS_FIT_TESTS=true
fi
"$RS" $W/dev/review-reslot-runtest.R $arm $path $pkg > $log 2>&1
grep RESULT $log || echo "NORESULT $arm $path"

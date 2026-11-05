#!/usr/bin/env bash
# Short check of fix round 2: the review's merged checks on the merged build.
cd /c/Users/adf44/source/r/frmtmb-wt-predread
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/review-predread-out/r2
export REV_LIB=C:/Users/adf44/source/r/rev-predread-lib/merge
"$RS" dev/review-predread-mergecheck.R > $O/mergecheck.log 2>&1
"$RS" dev/review-predread-r1-check4.R > $O/check4-merge.log 2>&1
"$RS" dev/review-predread-r1-check7.R > $O/check7-merge.log 2>&1
"$RS" dev/review-predread-r2-na.R > $O/na-merge.log 2>&1
"$RS" dev/review-predread-r2-cat.R > $O/cat-merge.log 2>&1
echo finished > $O/merge4.done

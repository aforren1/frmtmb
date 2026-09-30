#!/bin/bash
# Reviewer, punch round 1: every core and frmtmb.sample test file, gated,
# both arms, one file per process, 22 at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-formula2
O=dev/formula2-rev-p1
{
while read f; do for arm in before after; do echo "$arm frmtmb $f"; done; done < $O/core-files.txt
while read f; do for arm in before after; do echo "$arm frmtmb.sample $f"; done; done < dev/formula2-rev-tests/sample-files.txt
} | xargs -P 22 -n 3 sh -c \
  '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/formula2-rev-runtest-gated.R "$0" "$1" "$2" TRUE > "dev/formula2-rev-p1/$0-$1-$2.log" 2>&1'
echo ALLDONE

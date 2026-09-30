#!/bin/bash
# Reviewer driver: every frmtmb.sample test file, both arms, ungated.
cd /c/Users/adf44/source/r/frmtmb-wt-formula2
while read f; do
  for arm in before after; do echo "$arm $f"; done
done < dev/formula2-rev-tests/sample-files.txt | xargs -P 20 -n 2 sh -c \
  '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/formula2-rev-runtest-gated.R "$0" frmtmb.sample "$1" FALSE > "dev/formula2-rev-tests/s-$0-$1.log" 2>&1'
echo ALLDONE

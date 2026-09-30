#!/bin/bash
# Reviewer driver: every listed core test file, both arms, ungated, one
# file per process, 22 at a time.
cd /c/Users/adf44/source/r/frmtmb-wt-formula2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
while read f; do
  for arm in before after; do
    echo "$arm $f"
  done
done < dev/formula2-rev-tests/files.txt | xargs -P 22 -n 2 sh -c \
  '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/formula2-rev-runtest.R "$0" frmtmb "$1" FALSE > "dev/formula2-rev-tests/$0-$1.log" 2>&1'
echo ALLDONE

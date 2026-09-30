#!/usr/bin/env bash
# Punch round 1: every tier, in turn. Each script writes its own logs.
WT=/c/Users/adf44/source/r/frmtmb-wt-postfit2
cd "$WT"
rm -rf dev/postfit2-suite dev/postfit2-gated dev/postfit2-p1-sample dev/postfit2-p1-ext
bash dev/postfit2-suite.sh 2>&1 | grep -v "Segmentation fault" > dev/postfit2-log/p1-suite-summary.txt
bash dev/postfit2-gated.sh 2>&1 | grep -v "Segmentation fault" > dev/postfit2-log/p1-gated-tier.txt
bash dev/postfit2-sample-suite.sh 2>&1 | grep -v "Segmentation fault" > dev/postfit2-log/p1-sample-suite.txt
bash dev/postfit2-ext-suite.sh 2>&1 | grep -v "Segmentation fault" > dev/postfit2-log/p1-ext-suite.txt
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for f in tests/testthat/test-brms-suite-*.R; do
  b=$(basename "$f")
  "$R" dev/postfit2-runtest.R frmtmb "$b" --gated > dev/postfit2-log/p1-gated-core-$b.txt 2>&1 &
done
wait
echo ALLDONE > dev/postfit2-log/p1-all-done.txt

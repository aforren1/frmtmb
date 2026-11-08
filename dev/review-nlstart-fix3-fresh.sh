#!/bin/bash
# Review of lane nlstart, third short check: fresh designs on the
# fix build, both scripts, with the reference, two processes.
cd /c/Users/adf44/source/r/frmtmb-wt-nlstart
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
O=dev/review-nlstart-out
mkdir -p $O/fix3-fresh-after $O/fix3-fresh2-after
( for fs in gaussian,poisson,bernoulli binomial; do
    REV_ARM=after "$R" dev/review-nlstart-fresh.R $fs ref > "$O/fix3-fresh-after/$fs.log" 2>&1
  done ) &
( for fs in gamma,cumulative; do
    REV_ARM=after "$R" dev/review-nlstart-fresh.R $fs ref > "$O/fix3-fresh-after/$fs.log" 2>&1
  done
  REV_BODIES=logis3,gomp,sig4 REV_ARM=after "$R" dev/review-nlstart-fresh2.R gaussian,poisson,bernoulli,binomial,gamma,cumulative ref > "$O/fix3-fresh2-after/all.log" 2>&1 ) &
wait
echo ALLDONE

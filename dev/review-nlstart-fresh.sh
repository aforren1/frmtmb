#!/bin/bash
# Review of lane nlstart: fresh designs, one arm, three processes.
# Usage: dev/review-nlstart-fresh.sh <arm> <outdir> [ref]
arm=$1; out=$2; ref=$3
cd /c/Users/adf44/source/r/frmtmb-wt-nlstart
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
mkdir -p "$out"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for fs in gaussian,poisson bernoulli,binomial gamma,cumulative; do
  REV_ARM=$arm "$R" dev/review-nlstart-fresh.R $fs $ref > "$out/$fs.log" 2>&1 &
done
wait
echo ALLDONE

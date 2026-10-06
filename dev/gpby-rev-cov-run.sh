#!/bin/sh
# Reviewer coverage driver: one R process per (arm, mode, seed block).
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd /c/Users/adf44/source/r/frmtmb-wt-gpby
"$R" dev/gpby-rev-cov.R lane plugin 1 100 &
"$R" dev/gpby-rev-cov.R base plugin 1 100 &
"$R" dev/gpby-rev-cov.R lane oracle 1 100 &
"$R" dev/gpby-rev-cov.R lane oracle 101 200 &
"$R" dev/gpby-rev-cov.R base oracle 1 100 &
"$R" dev/gpby-rev-cov.R base oracle 101 200 &
"$R" dev/gpby-crit.R base > dev/gpby-rev-crit-base.txt 2>&1 &
"$R" dev/gpby-crit.R lane > dev/gpby-rev-crit-lane.txt 2>&1 &
wait
echo ALLDONE

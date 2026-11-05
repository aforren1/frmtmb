#!/usr/bin/env bash
# Final runs, first track: the grid, the review's off-grid designs,
# the year and scale tables, then the suites one file per process.
cd /c/Users/adf44/source/r/frmtmb-wt-impsamp || exit 1
L=dev/impsamp-out
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LIB=C:/Users/adf44/source/r/wt-impsamp-lib
bash dev/impsamp-grid-run.sh final "$LIB" > $L/grid-final-run.txt 2>&1
"$R" dev/impsamp-grid-summary.R base final > $L/grid-summary-final.txt 2>&1
echo grid-done > $L/final-a-grid.txt
bash dev/impsamp-offgrid.sh 6 > $L/offgrid-run.txt 2>&1
"$R" dev/impsamp-offgrid-sum.R > $L/offgrid-summary.txt 2>&1
"$R" dev/impsamp-year.R after > $L/year-after.log 2>&1
"$R" dev/review-impsamp-scale-se.R > $L/scale-se.txt 2>&1
"$R" dev/impsamp-scales.R > $L/scales.txt 2>&1
"$R" dev/impsamp-natural.R > $L/natural.txt 2>&1
"$R" dev/impsamp-classes.R > $L/classes.txt 2>&1
"$R" dev/impsamp-ordfix.R > $L/ordfix.txt 2>&1
"$R" dev/impsamp-escape.R > $L/escape-after.txt 2>&1
echo tables-done > $L/final-a-tables.txt
rm -rf $L/suite-frmtmb $L/suite-frmtmb.learn $L/suite-frmtmb.spatial $L/suite-frmtmb.sample $L/suite-frmtmb.latent $L/suite-frmtmb.eam
bash dev/impsamp-suite.sh frmtmb tests/testthat 6 > $L/suite-core.txt 2>&1
for p in learn spatial sample latent eam; do
  bash dev/impsamp-suite.sh frmtmb.$p extensions/frmtmb.$p/tests/testthat 6 > $L/suite-$p.txt 2>&1
done
echo suites-done > $L/final-a-suites.txt

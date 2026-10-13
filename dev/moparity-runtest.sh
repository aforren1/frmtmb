#!/bin/sh
# Run one test file in its own R process against a library set.
# Usage: sh dev/moparity-runtest.sh <pkg> <file> <log> [base]
# base: run against rellib-r7 only (the "before" arm)
pkg=$1; f=$2; log=$3; arm=$4
WT=C:/Users/adf44/source/r/frmtmb-wt-moparity
if [ "$pkg" = "frmtmb" ]; then dir=$WT/tests/testthat; else dir=$WT/extensions/$pkg/tests/testthat; fi
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export MOPARITY_ARM=$arm
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" -e "
arm <- Sys.getenv('MOPARITY_ARM')
libs <- c('C:/Users/adf44/source/r/wt-moparity-lib', 'C:/Users/adf44/source/r/rellib-r7', 'C:/Users/adf44/AppData/Local/R/win-library/4.6')
if (identical(arm, 'base')) libs <- libs[-1]
.libPaths(libs)
p <- '$pkg'
suppressPackageStartupMessages(library(p, character.only = TRUE))
cat('lib:', find.package('frmtmb'), find.package(p), '\n')
setwd('$dir')
r <- testthat::test_file('$f', package = p, env = testthat::test_env(p), reporter = 'summary', stop_on_failure = FALSE)
d <- as.data.frame(r)
cat(sprintf('RESULT %s pass=%d fail=%d err=%d skip=%d warn=%d\n', '$f', sum(d\$passed), sum(d\$failed), sum(d\$error), sum(d\$skipped), sum(d\$warning)))
" > "$log" 2>&1
grep "RESULT\|^lib:" "$log"

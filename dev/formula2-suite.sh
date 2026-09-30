#!/usr/bin/env bash
# Run every test file of a package, one R process per file, P at a
# time. Args: <pkg> <mode: plain|gated> <P>. Logs in dev/formula2-suite/.
pkg=$1; mode=$2; P=${3:-12}
wt=/c/Users/adf44/source/r/frmtmb-wt-formula2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
if [ "$pkg" = frmtmb ]; then dir=$wt/tests/testthat; else dir=$wt/extensions/$pkg/tests/testthat; fi
runner=$wt/dev/formula2-run-test.R
[ "$mode" = gated ] && runner=$wt/dev/formula2-run-gated.R
ls $dir | grep '^test-.*\.R$' | xargs -P $P -I{} sh -c "\"$R\" $runner $pkg {} > $wt/dev/formula2-suite/$pkg-$mode-{}.log 2>&1"

#!/bin/sh
# Lane sampfix: every test file of one package, one R process per file,
# P at a time.   sh dev/sampfix-suite.sh <lane|ref> <pkg> <outdir> <P>
arm=$1; pkg=$2; out=$3; P=$4
WT=/c/Users/adf44/source/r/frmtmb-wt-sampfix
if [ "$pkg" = "frmtmb" ]; then dir=$WT/tests/testthat; else dir=$WT/extensions/$pkg/tests/testthat; fi
mkdir -p "$out"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/stan-cache
ls "$dir" | grep '^test-.*\.R$' | xargs -P "$P" -I{} sh -c \
  '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" '"$WT"'/dev/sampfix-runtest.R '"$arm $pkg"' {} > '"$out"'/{}.txt 2>&1'
grep -h "^RESULT" "$out"/*.txt

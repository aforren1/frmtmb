#!/bin/sh
# Reviewer, lane sampfix: run test files one R process each, P at a time.
#   sh dev/sampfix-rev-suite.sh <lane|ref> <pkg> <outdir> <P> [file ...]
arm=$1; pkg=$2; out=$3; P=$4; shift 4
WT=/c/Users/adf44/source/r/frmtmb-wt-sampfix
if [ "$pkg" = "frmtmb" ]; then dir=$WT/tests/testthat; else dir=$WT/extensions/$pkg/tests/testthat; fi
mkdir -p "$out"
if [ $# -gt 0 ]; then files="$*"; else files=$(ls "$dir" | grep '^test-.*\.R$'); fi
printf '%s\n' $files | xargs -P "$P" -I{} sh -c \
  '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" '"$WT"'/dev/sampfix-rev-runtest.R '"$arm $pkg"' {} > '"$out"'/{}.txt 2>&1'
grep -h "^RESULT" "$out"/*.txt

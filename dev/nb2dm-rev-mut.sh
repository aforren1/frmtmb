#!/usr/bin/env bash
# Reviewer: each mutation of dev/nb2dm-rev-mut.R against the test files
# meant to catch it, one R process per (mutation, file), two at a time.
cd "$(dirname "$0")/.."
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=dev/nb2dm-rev-out/mut
mkdir -p $OUT
T=tests/testthat
S=extensions/frmtmb.sample/tests/testthat
jobs_list="none:$T/test-negbinomial2.R none:$T/test-dirichlet-multinomial.R
mnfam:$T/test-dirichlet-multinomial.R mnfam:$S/test-nb2dm-draws.R
mnfam:$T/test-simulate-density.R
noeta:$T/test-negbinomial2.R nopost:$T/test-negbinomial2.R
plain:$T/test-dirichlet-multinomial.R nostart:$T/test-negbinomial2.R"
for j in $jobs_list; do
  m=${j%%:*}; f=${j#*:}
  while [ "$(jobs -rp | wc -l)" -ge 2 ]; do sleep 1; done
  ( "$RS" dev/nb2dm-rev-mut.R "$m" "$f" > "$OUT/$m--$(basename "$f" .R).txt" 2>&1 ) &
done
wait
grep -h -A12 "^RESULT" $OUT/*.txt > dev/nb2dm-rev-out/mut.txt

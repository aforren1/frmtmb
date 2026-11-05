#!/bin/bash
# Runs every hand-mapped identity case, 8 at once.
cd /c/Users/adf44/source/r/frmtmb-wt-latres
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for c in ma arma cosy unstr ar_uns ma_uns unstr_uns ar_rag arma_rag cosy_rag ar_gap ar_single ar_notime ar_nogr cosy_nogr ma_notime ar_gamma ar_bern condar condar2 condar_gamma condar_gamma_inv condar_probit; do
  echo $c
done | xargs -P 8 -I{} sh -c "\"$R\" dev/review-latres-id.R {} dev/review-latres-out/id-{}.tsv > dev/review-latres-out/id-{}.log 2>&1"

#!/usr/bin/env bash
# M1 and the warning on the warm-started build: the reviewer's floor and noise-free
# scripts and gpgap's reviewer's floor script, read-only from that tree.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
O=dev/ncgp-out/m1
mkdir -p $O
"$RS" dev/review-ncgp-floor.R lane 40 > $O/floor-40.txt 2>&1 &
"$RS" dev/review-ncgp-floor.R lane 100 > $O/floor-100.txt 2>&1 &
"$RS" dev/review-ncgp-nfpred.R > $O/nfpred.txt 2>&1 &
"$RS" dev/review-ncgp-mpfr.R > $O/mpfr.txt 2>&1 &
"$RS" dev/review-ncgp-lanerows.R > $O/lanerows.txt 2>&1 &
"$RS" dev/review-ncgp-stationary.R > $O/stationary.txt 2>&1 &
"$RS" ../frmtmb-wt-gpgap/dev/review-gpgap-ncgpfloor.R > $O/ncgpfloor.txt 2>&1 &
wait
echo M1-DONE

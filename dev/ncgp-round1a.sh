#!/usr/bin/env bash
# Fix round 1's measurements, on the final build. Bounded at 8 R processes.
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r/wt-ncgp-lib
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
O=dev/ncgp-out
mkdir -p $O/m1 $O/rs
"$RS" dev/review-ncgp-modes2.R lane > $O/modes2-lane.txt 2>&1 &
"$RS" dev/ncgp-equiv.R lane > $O/equiv-lane.txt 2>&1 &
"$RS" dev/ncgp-keep.R > $O/keep.txt 2>&1 &
"$RS" dev/ncgp-truth.R lane $L > $O/truth-lane.txt 2>&1 &
"$RS" dev/ncgp-bigunit.R base > $O/bigunit-base.txt 2>&1 &
"$RS" dev/ncgp-bigunit.R $L > $O/bigunit-lane.txt 2>&1 &
wait
"$RS" dev/ncgp-nuggetrate.R > $O/nuggetrate.txt 2>&1 &
"$RS" dev/review-ncgp-floor.R lane 40 > $O/m1/floor-40.txt 2>&1 &
"$RS" dev/review-ncgp-floor.R lane 100 > $O/m1/floor-100.txt 2>&1 &
"$RS" dev/review-ncgp-nfpred.R > $O/m1/nfpred.txt 2>&1 &
"$RS" dev/review-ncgp-mpfr.R > $O/m1/mpfr.txt 2>&1 &
"$RS" dev/review-ncgp-lanerows.R > $O/m1/lanerows.txt 2>&1 &
"$RS" dev/review-ncgp-stationary.R > $O/m1/stationary.txt 2>&1 &
"$RS" ../frmtmb-wt-gpgap/dev/review-gpgap-ncgpfloor.R > $O/m1/ncgpfloor.txt 2>&1 &
wait
"$RS" dev/ncgp-nuggetrate-sum.R > $O/nuggetrate-sum.txt 2>&1
"$RS" dev/ncgp-randstart.R 1e-6 default > $O/rs/1e-6-default.txt 2>&1 &
"$RS" dev/ncgp-randstart.R 1e-6 mode > $O/rs/1e-6-mode.txt 2>&1 &
"$RS" dev/ncgp-randstart.R 1e-3 default > $O/rs/1e-3-default.txt 2>&1 &
wait
"$RS" dev/review-ncgp-mem.R lane 300 2 > $O/mem-lane-300.txt 2>&1
echo ROUND1A-DONE

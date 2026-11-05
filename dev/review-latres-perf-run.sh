#!/bin/bash
cd /c/Users/adf44/source/r/frmtmb-wt-latres
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/review-latres-out/perf
rm -f dev/review-latres-out/perf/*
run() { "$R" dev/review-latres-perf.R "$1" "$2" "$3" dev/review-latres-out/perf/$4.tsv > dev/review-latres-out/perf/$4.log 2>&1; }
run 1 2000 "ar(t, g, cov = TRUE)" a1
run 1 2000 "ma(t, g, cov = TRUE)" a2
run 1 2000 "arma(t, g, cov = TRUE)" a3
run 1 2000 "ar(t, g, p = 2)" a4
run 1 2000 "cosy(t, g)" a5
run 20 100 "ar(t, g, cov = TRUE)" a6
run 400 20 "ar(t, g, cov = TRUE)" a7
run 500 10 "unstr(t, g)" u1
run 500 20 "unstr(t, g)" u2
run 20 100 "ar1(ft + 0 | g)" d1
run 50 20 "toep(ft + 0 | g)" d2
run 1 20000 "ar(t, g, cov = TRUE)" b1
run 1 20000 "arma(t, g, cov = TRUE)" b2
run 1 20000 "ar(t, g, p = 2)" b3
run 5000 4 "ar(t, g, cov = TRUE)" c1
run 200 100 "ar(t, g, cov = TRUE)" c2
run 1 100000 "ar(t, g, cov = TRUE)" e1

#!/usr/bin/env bash
# The calibration of the warning on the review's two tables.
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin"
mkdir -p dev/refuse-f2-out
for f in Gamma exponential beta exgaussian student lognormal von_mises gaussian; do echo $f; done |
  xargs -P 4 -I{} bash -c '"$0" dev/refuse-f2-laplace.R {} > dev/refuse-f2-out/lap-{}.log 2>&1' "$R"

#!/usr/bin/env bash
f="$1"
b=$(basename "$f")
R="C:/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=dev/resmooth-rev-suite
case "$b" in
  test-brms-likelihood.R|test-brms-methods.R) export REV_GATE=off ;;
esac
"$R" dev/resmooth-rev-run-one.R "$f" > "$OUT/$b.txt" 2>&1
echo "done $b gate=${REV_GATE:-on}" >> "$OUT/driver.log"

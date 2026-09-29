#!/usr/bin/env bash
f="$1"; b=$(basename "$f")
"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/resmooth-rev-run-one.R "$f" > "dev/resmooth-rev2-suite/$b.txt" 2>&1
echo "done $b" >> dev/resmooth-rev2-suite/driver.log

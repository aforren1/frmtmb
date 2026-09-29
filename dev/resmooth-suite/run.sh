#!/bin/sh
cd /c/Users/adf44/source/r/frmtmb-wt-resmooth
cat dev/resmooth-suite-files.lst | xargs -P 8 -I{} sh -c '
  b=$(basename "{}")
  "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/resmooth-run-one.R "{}" > "dev/resmooth-suite/$b.txt" 2>&1
'
echo SUITE-DONE

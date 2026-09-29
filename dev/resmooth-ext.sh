#!/bin/sh
cd /c/Users/adf44/source/r/frmtmb-wt-resmooth
ls extensions/frmtmb.sample/tests/testthat/test-*.R \
   extensions/frmtmb.spline/tests/testthat/test-*.R \
   > dev/resmooth-extfiles.lst
mkdir -p dev/resmooth-ext
cat dev/resmooth-extfiles.lst | xargs -P 3 -I{} sh -c '
  b=$(echo "{}" | sed "s|extensions/||; s|/tests/testthat/|-|")
  "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/resmooth-run-one.R "{}" > "dev/resmooth-ext/$b.txt" 2>&1
'
echo EXT-DONE

#!/bin/sh
cd /c/Users/adf44/source/r/frmtmb-wt-resmooth
cat dev/resmooth-extfiles.lst | xargs -P 6 -I{} sh -c '
  b=$(echo "{}" | sed "s|extensions/||; s|/tests/testthat/|-|")
  "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/resmooth-run-one.R "{}" > "dev/resmooth-ext2/$b.txt" 2>&1
'
echo EXT2-DONE

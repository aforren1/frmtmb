#!/usr/bin/env bash
# Lane drawdir, fix round 1 with option (b): the 12 random-slope fires
# again (dev/drawdir-fires12.R), then the draw test files under
# OpenBLAS 0.3.26 with LAPACK (dev/drawdir-ob2.sh). One R process.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-drawdir
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/drawdir-fires12.R slpw:12 slpw:14 slpw:15 slpw:18 slp10:10 prsp:1 prsp:10 ordrs:9 S0.1:14 S0.1:16 S0.1:24 S0.2:8 > dev/drawdir-fires12b.txt 2>&1
rm -rf dev/drawdir-ob2-out
bash dev/drawdir-ob2.sh > dev/drawdir-ob2.txt 2>&1
rm -rf dev/ciharden-out
echo batchC done

#!/usr/bin/env bash
# Lane drawdir: the draw test files on rellib-r22 (seen to fail), then
# every suite (core and the eight extensions), one file per process,
# one process at a time (the other slot runs dev/drawdir-batchA.sh).
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-drawdir
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true R_PARALLEL_PORT=11843
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p dev/drawdir-base-out
for f in test-predict-draw-check test-predict-draw-edge; do
  DRAWDIR_BASE=1 "$RS" dev/drawdir-run-tests.R frmtmb tests/testthat/$f.R \
    > dev/drawdir-base-out/$f.txt 2>&1
done
bash dev/drawdir-run-par.sh all dev/drawdir-jobs-all.lst 1
echo batchB done

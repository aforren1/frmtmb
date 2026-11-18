#!/usr/bin/env bash
# Lane drawdir: the measurements on the lane build, one R process at a
# time: the verdict on every group (dev/drawdir-verdict.R), the short
# check's script on the lane library (dev/drawdir-repro.R), the ulp runs
# of the new tests' fits (dev/drawdir-ulp.R), and the cost
# (dev/drawdir-cost.R).
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-drawdir
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LIB=C:/Users/adf44/source/r/wt-drawdir-lib
for g in gglm gglmm wb kd zip sweep target mix40 partR; do
  "$RS" dev/drawdir-verdict.R $LIB dev/drawdir-v-$g.tsv $g > dev/drawdir-v-$g.txt 2>&1
done
DRAWDIR_LIB=$LIB "$RS" dev/drawdir-repro.R > dev/drawdir-repro-lane.txt 2>&1
"$RS" dev/drawdir-ulp.R $LIB ref > dev/drawdir-ulp-ref.txt 2>&1
"$RS" dev/drawdir-cost.R $LIB > dev/drawdir-cost-lane.txt 2>&1
echo batchA done

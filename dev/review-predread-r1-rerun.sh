#!/usr/bin/env bash
# Short check of fix round 1: the review's own scripts on the lane's final build,
# logs under dev/review-predread-out/r1/. N at a time.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-predread
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
O=dev/review-predread-out/r1
mkdir -p "$O"
N=${1:-4}
run() {
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done
  ( "$RS" "dev/review-predread-$1.R" ${2:-} > "$O/$1.log" 2>&1 ) &
}
run mefam
run me
run me-edge
run misc
run resp
run cells
run mejoint
run share
run update
run update2
run regress dev/review-predread-out/regress-after.rds
run regress-me dev/review-predread-out/regme-after.rds
run explain
wait
"$RS" dev/review-predread-regress-cmp.R > "$O/regress-cmp.log" 2>&1
"$RS" dev/review-predread-out/regme-cmp.R > "$O/regme-cmp.log" 2>&1
echo finished > "$O/rerun.done"

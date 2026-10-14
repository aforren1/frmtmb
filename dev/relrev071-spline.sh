#!/usr/bin/env bash
# Release reviewer of 0.71.0: frmtmb.spline's floor in both directions.
# new: frmtmb.spline 0.10.2 (this tree) over frmtmb 0.70.0 (rellib-r8),
#      the tree's spline tests;
# old: frmtmb.spline 0.10.1 (rellib-r8's) over frmtmb 0.71.0 (rellib-r9),
#      HEAD's spline tests plus the tree's new test-illcond.R.
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
W=$REL/dev/relrev071-log/spline
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
L=C:/Users/adf44/source/r
export NOT_CRAN=true
rm -rf "$W"; mkdir -p "$W"
N=8
jobs_wait() { while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 2; done; }
run() { jobs_wait; "$RS" "$REL/dev/relrev071-runtest.R" "$@" < /dev/null \
  > "$W/$TAG-$(basename "$3" .R).txt" 2>&1 & }
TAG=new
for f in "$REL"/extensions/frmtmb.spline/tests/testthat/test-*.R; do
  run "$L/relrev071-lib/newspline;$L/rellib-r8" frmtmb.spline "$f"
done
TAG=old
for f in "$REL"/dev/relrev071-out/oldspline-tests/test-*.R \
         "$REL"/extensions/frmtmb.spline/tests/testthat/test-illcond.R; do
  run "$L/relrev071-lib/oldspline;$L/rellib-r9" frmtmb.spline "$f"
done
wait
grep -h RESULT "$W"/*.txt > /dev/null
for t in new old; do
  echo "== $t"; grep -H '^RESULT' "$W"/$t-*.txt | sed 's|.*/||'
done > "$W/summary.txt"

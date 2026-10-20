#!/usr/bin/env bash
# Reviewer, check of fix round 2: the test file, the new probe, and the
# earlier probes rerun, one R process each, on the lane build.
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drop1raw
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/drop1raw-rev-out"
export NOT_CRAN=true
cd "$ROOT"
"$RS" dev/drop1raw-runtest.R lane frmtmb tests/testthat/test-refit-data.R > "$OUT/r2-t-refit-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-r2.R > "$OUT/r2-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-r1.R lane > "$OUT/r2-r1-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-r1e.R > "$OUT/r2-r1e-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-probe.R lane > "$OUT/r2-probe-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-step.R lane > "$OUT/r2-step-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-step2.R lane > "$OUT/r2-step2-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-locals.R lane > "$OUT/r2-locals-lane.txt" 2>&1 &
wait
echo done

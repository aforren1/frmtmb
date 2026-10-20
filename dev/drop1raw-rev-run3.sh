#!/usr/bin/env bash
# Reviewer, re-check of fix round 1: the earlier probes, the new probe
# and the two pinning test files, one R process each, on the lane build.
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drop1raw
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/drop1raw-rev-out"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
cd "$ROOT"
"$RS" dev/drop1raw-runtest.R lane frmtmb tests/testthat/test-refit-data.R \
  > "$OUT/r1-t-refit-lane.txt" 2>&1 &
"$RS" dev/drop1raw-runtest.R lane frmtmb.latent \
  extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R \
  > "$OUT/r1-t-hmm-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-r1.R lane > "$OUT/r1-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-probe.R lane > "$OUT/r1-probe-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-step.R lane > "$OUT/r1-step-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-step2.R lane > "$OUT/r1-step2-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-locals.R lane > "$OUT/r1-locals-lane.txt" 2>&1 &
"$RS" dev/drop1raw-rev-hmmstart.R lane > "$OUT/r1-hmmstart-lane.txt" 2>&1 &
wait
echo done

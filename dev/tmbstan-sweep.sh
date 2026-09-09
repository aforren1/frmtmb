#!/bin/sh
# lane tmbstan: run every touched test file, one file per R process,
# in the arm named by $1 ("clean" or "broken"). Every count field is
# printed, errors included: a runner that sums `failed` alone reports a
# clean line for a file that aborted, which is exactly what the broken
# arm does.
ARM="${1:-clean}"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
S="C:/Users/adf44/source/r/frmtmb-wt-tmbstan/dev/tmbstan-run.R"
FILES="test-tmbstan-build-guard.R test-conditional-effects-draws.R \
test-draws-methods.R test-draws-spellings.R test-evidence-ratio.R \
test-loo.R test-parallel-chains.R test-reparam.R test-sample-direct.R \
test-sampling-ported.R test-scale.R test-stan-control.R \
test-compat-preflight.R"
for f in $FILES; do
  if [ "$ARM" = "broken" ]; then
    "$R" "$S" "$f" "" broken 2>&1 | grep '^###'
  else
    "$R" "$S" "$f" 2>&1 | grep '^###'
  fi
done

#!/bin/sh
# Reviewer: run the listed test files, one per process, 4 at a time.
W=/c/Users/adf44/source/r/frmtmb-wt-priordraws
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LOG=$W/dev/review-priordraws-log
export FRMTMB_BRMS_FIT_TESTS=true
run() {
  # $1 tag, $2 pkg, $3 file, $4 arm
  "$RS" $W/dev/review-priordraws-runtest.R "$2" "$3" "$4" > "$LOG/t-$1.txt" 2>&1
  echo "done $1"
}
C=$W/tests/testthat
S=$W/extensions/frmtmb.sample/tests/testthat
run core-prioronly frmtmb $C/test-prior-only-objective.R lane &
run core-msguniq frmtmb $C/test-message-uniqueness.R lane &
run core-bracket frmtmb $C/test-bracket-access.R lane &
run s-priordraws1 frmtmb.sample $S/test-prior-draws.R lane &
wait
run s-priordraws2 frmtmb.sample $S/test-prior-draws.R lane &
run s-priordraws3 frmtmb.sample $S/test-prior-draws.R lane &
run s-warmup frmtmb.sample $S/test-warmup-draws.R lane &
run s-evid frmtmb.sample $S/test-evidence-ratio.R lane &
wait
run s-methods frmtmb.sample $S/test-draws-methods.R lane &
run s-spell frmtmb.sample $S/test-draws-spellings.R lane &
run s-loo frmtmb.sample $S/test-loo.R lane &
run s-generic frmtmb.sample $S/test-generic-collision.R lane &
wait
run s-brmssuite frmtmb.sample $S/test-brms-suite-methods.R lane &
run s-bracket frmtmb.sample $S/test-bracket-access.R lane &
run s-msguniq frmtmb.sample $S/test-message-uniqueness.R lane &
run core-prioronly-base frmtmb $C/test-prior-only-objective.R base &
wait
echo ALLDONE

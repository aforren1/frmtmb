#!/bin/sh
# Reviewer, short check of fix round 1: one file per process, 4 at a time.
W=/c/Users/adf44/source/r/frmtmb-wt-priordraws
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
LOG=$W/dev/review-priordraws-log
export FRMTMB_BRMS_FIT_TESTS=true
run() { "$RS" $W/dev/review-priordraws-runtest.R "$2" "$3" lane > "$LOG/fr1-$1.txt" 2>&1; }
C=$W/tests/testthat
S=$W/extensions/frmtmb.sample/tests/testthat
run pd1 frmtmb.sample $S/test-prior-draws.R & run pd2 frmtmb.sample $S/test-prior-draws.R &
run pd3 frmtmb.sample $S/test-prior-draws.R & run er1 frmtmb.sample $S/test-evidence-ratio.R &
wait
run er2 frmtmb.sample $S/test-evidence-ratio.R & run er3 frmtmb.sample $S/test-evidence-ratio.R &
run wu frmtmb.sample $S/test-warmup-draws.R & run dm frmtmb.sample $S/test-draws-methods.R &
wait
run cpo frmtmb $C/test-prior-only-objective.R & run cme frmtmb $C/test-me.R &
run cmu frmtmb $C/test-message-uniqueness.R & run smu frmtmb.sample $S/test-message-uniqueness.R &
wait
echo ALLDONE

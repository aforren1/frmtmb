#!/usr/bin/env bash
# Lane gpselimit, fix round 1: every suite the fixes reach, in order, two R
# processes at most: all 550 test files ungated, the 45 gp() and SE files
# gated, the one-session runs (core forward; frmtmb.sample three ways),
# then R CMD check of core and of frmtmb.sample.
# Usage: bash dev/gpselimit-fx-chain.sh
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-gpselimit
O=dev/gpselimit-out
bash dev/gpselimit-par.sh fx-all $O/jobs-all.txt 2 > /dev/null 2>&1
FRMTMB_BRMS_FIT_TESTS=true bash dev/gpselimit-par.sh fx-t1 $O/jobs-t1.txt 2 \
  > /dev/null 2>&1
rm -f $O/sessions.log
bash dev/gpselimit-sessions.sh frmtmb:forward frmtmb.sample:forward \
  frmtmb.sample:reverse frmtmb.sample:preload > /dev/null 2>&1
bash dev/gpselimit-check.sh > $O/check-summary-fx.txt 2>&1
cp dev/gpselimit-check/frmtmb.Rcheck/00check.log $O/00check-frmtmb-fx.log
bash dev/gpselimit-check-sample.sh > $O/check-sample-summary-fx.txt 2>&1
cp dev/gpselimit-check/frmtmb.sample/frmtmb.sample.Rcheck/00check.log \
  $O/00check-frmtmb.sample-fx.log
cp dev/gpselimit-check/frmtmb.Rcheck/tests/testthat.Rout.fail \
  $O/ 2>/dev/null
grep -a "\[ FAIL" dev/gpselimit-check/frmtmb.Rcheck/tests/testthat.Rout \
  dev/gpselimit-check/frmtmb.sample/frmtmb.sample.Rcheck/tests/testthat.Rout \
  > $O/check-tests-fx.txt 2>&1
echo "CHAIN DONE" > $O/fx-chain.done

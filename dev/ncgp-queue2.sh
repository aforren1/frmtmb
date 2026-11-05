#!/usr/bin/env bash
# After the gated tier: one-session runs of all nine suites, ungated, on this
# build, then OpenBLAS, then R CMD check of the two changed packages.
set -u
until grep -q SEEN-DONE dev/ncgp-out/queue.log; do sleep 20; done
export NOT_CRAN=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-ncgp-lib
export TMP="${TMP:-$TEMP}"
bash dev/release/run-sessions.sh ncgpsess 6 > dev/ncgp-out/sess-driver.txt 2>&1
cp dev/release/ncgpsess.log dev/ncgp-out/ncgpsess.log 2>/dev/null
echo SESS-DONE
bash dev/ncgp-ob.sh > dev/ncgp-out/ob.log 2>&1
echo OB-DONE
bash dev/ncgp-check.sh . > dev/ncgp-out/check-core.txt 2>&1
bash dev/ncgp-check.sh extensions/frmtmb.sample > dev/ncgp-out/check-sample.txt 2>&1
echo CHECK-DONE

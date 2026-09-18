#!/bin/sh
# The lane's final pass: the core suite, all seven extension suites, the
# gated tiers and the ported brms tier, one R process per file, nothing
# else running.
#
#   sh dev/adefects-all.sh
#
# It writes dev/adefects-log/runs-complete.txt last, so a waiter can
# poll for that file rather than for a process.
set -e
cd /c/Users/adf44/source/r/frmtmb-wt-adefects
rm -f dev/adefects-log/runs-complete.txt
sh dev/adefects-suite.sh frmtmb tests/testthat \
  > dev/adefects-log/suite-core.txt 2>&1
for p in coupling eam latent learn ode sample spline; do
  sh dev/adefects-suite.sh "frmtmb.$p" "extensions/frmtmb.$p/tests/testthat" \
    > "dev/adefects-log/suite-$p.txt" 2>&1
done
sh dev/adefects-gated.sh > dev/adefects-log/gated.txt 2>&1
FRMTMB_PORT_ROOT=/c/Users/adf44/source/r/frmtmb-wt-adefects \
  FRMTMB_PORT_LIB=C:/Users/adf44/source/r/adefects-lib \
  sh dev/brmsport-tier.sh > dev/adefects-log/port-tier.txt 2>&1
echo ALLDONE > dev/adefects-log/runs-complete.txt

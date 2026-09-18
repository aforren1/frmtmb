#!/bin/sh
# The gated tiers, the job list of dev/release/run-gated.ps1, against
# this lane's library. One R process per file.
#
#   sh dev/adefects-gated.sh > dev/adefects-log/gated.txt 2>&1
#
# The brms-suite files are run by dev/brmsport-tier.sh instead, which
# asserts their verdicts; they are left out here so the two runs do not
# double count.
set -e
cd /c/Users/adf44/source/r/frmtmb-wt-adefects
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_FUZZ=true
cache=/c/Users/adf44/source/r/frmtmb-wt-adefects/dev/stan-cache
export FRMTMB_STAN_CACHE="$cache"
rt=/c/rtools45
export PATH="$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin:$PATH"
core=tests/testthat
n=0
for f in "$core"/test-bcm-*.R \
         "$core/test-brms-agreement.R" "$core/test-brms-likelihood.R" \
         "$core/test-brms-methods.R" "$core/test-brms-priors.R" \
         "$core/test-brms-port.R" "$core/test-rl-example.R" \
         "$core/test-fuzz.R"; do
  n=$((n + 1))
  "$R" dev/adefects-run.R frmtmb "$f" 2>&1 |
    grep -E "^RESULT|^  BLOCK" || echo "RESULT frmtmb $f NOOUTPUT"
done
sdir=extensions/frmtmb.sample/tests/testthat
for f in "$sdir/test-loo.R" "$sdir/test-sampling-ported.R"; do
  if [ ! -f "$f" ]; then echo "MISSING $f"; exit 1; fi
  n=$((n + 1))
  "$R" dev/adefects-run.R frmtmb.sample "$f" 2>&1 |
    grep -E "^RESULT|^  BLOCK" || echo "RESULT frmtmb.sample $f NOOUTPUT"
done
lf=extensions/frmtmb.learn/tests/testthat/test-stan-identity.R
if [ ! -f "$lf" ]; then echo "MISSING $lf"; exit 1; fi
n=$((n + 1))
"$R" dev/adefects-run.R frmtmb.learn "$lf" 2>&1 |
  grep -E "^RESULT|^  BLOCK" || echo "RESULT frmtmb.learn $lf NOOUTPUT"
echo "GATED ran $n files"

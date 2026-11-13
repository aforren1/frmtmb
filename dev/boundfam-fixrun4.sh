#!/usr/bin/env bash
# Lane boundfam consolidation fixes: test-se-limit.R on rellib-r20, the
# core suite with frmtmb.eam and frmtmb.learn (traced), the gated core
# test-brms-*.R files, then R CMD check --as-cran of core.
set -u
cd "$(dirname "$0")/.."
export NOT_CRAN=true
bash dev/boundfam-par2.sh base dev/boundfam-jobs-sl.txt 1
rm -rf dev/boundfam-suite5-base; mv dev/boundfam-suite2-base dev/boundfam-suite5-base
bash dev/boundfam-par2.sh lane dev/boundfam-jobs-cons.txt 3
rm -rf dev/boundfam-suite5-lane; mv dev/boundfam-suite2-lane dev/boundfam-suite5-lane
FRMTMB_BRMS_FIT_TESTS=true bash dev/boundfam-par2.sh lane dev/boundfam-jobs-gated.txt 3
rm -rf dev/boundfam-suite5-gated; mv dev/boundfam-suite2-lane dev/boundfam-suite5-gated
unset NOT_CRAN
bash dev/boundfam-check.sh . > dev/boundfam-out/check-core-cons.txt 2>&1
echo "fixrun4 done"

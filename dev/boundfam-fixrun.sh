#!/usr/bin/env bash
# Lane boundfam fix round: the ungated run (core, frmtmb.eam,
# frmtmb.learn, frmtmb.spatial test-sar.R), traced, then the gated core
# test-brms-*.R files, each one file per process, three at a time.
set -u
cd "$(dirname "$0")/.."
export NOT_CRAN=true
bash dev/boundfam-par2.sh lane dev/boundfam-jobs-fix.txt 3
mv dev/boundfam-suite2-lane dev/boundfam-suite3-lane
FRMTMB_BRMS_FIT_TESTS=true bash dev/boundfam-par2.sh lane dev/boundfam-jobs-gated.txt 3
mv dev/boundfam-suite2-lane dev/boundfam-suite3-gated
echo "fixrun done"

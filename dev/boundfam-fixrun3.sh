#!/usr/bin/env bash
# Lane boundfam fix round 2 (held limit directions): the core suite,
# one file per process, traced; the gated core test-brms-*.R files; then
# R CMD check --as-cran of core.
set -u
cd "$(dirname "$0")/.."
export NOT_CRAN=true
bash dev/boundfam-par2.sh lane dev/boundfam-jobs-core318.txt 3
rm -rf dev/boundfam-suite4-lane; mv dev/boundfam-suite2-lane dev/boundfam-suite4-lane
FRMTMB_BRMS_FIT_TESTS=true bash dev/boundfam-par2.sh lane dev/boundfam-jobs-gated.txt 3
rm -rf dev/boundfam-suite4-gated; mv dev/boundfam-suite2-lane dev/boundfam-suite4-gated
unset NOT_CRAN
bash dev/boundfam-check.sh . > dev/boundfam-out/check-core-fix2.txt 2>&1
echo "fixrun3 done"

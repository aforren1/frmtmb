#!/usr/bin/env bash
# Lane silent2: every test file of core and the eight extensions, one R
# process per file, three at a time, on the lane library (then rellib-r21),
# with the gated brms tests on.
#   bash dev/silent2-full.sh <tier> <job list>
set -u
: "${TMP:?TMP is unset}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-silent2-lib
# R_PARALLEL_PORT stays unset: each process then draws its own port, so
# the three at once do not share one
unset R_PARALLEL_PORT
bash "$ROOT/dev/release/run-par.sh" "$1" "$2" 3

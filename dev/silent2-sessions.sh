#!/usr/bin/env bash
# Lane silent2: every suite in one session, three orders, on the lane
# library (then rellib-r21), N sessions at a time.
#   bash dev/silent2-sessions.sh <tier> <N>
set -u
: "${TMP:?TMP is unset}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-silent2-lib
unset R_PARALLEL_PORT
bash "$ROOT/dev/release/run-sessions.sh" "$1" "$2"

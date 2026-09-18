#!/bin/sh
# Run a package's whole test suite, ONE R PROCESS PER FILE. A single
# process has repeatedly hidden state leakage in this project, and a
# runner that sums `failed` and not `error` prints a clean line for a
# file that aborted halfway, so both are summed and printed.
#
#   sh dev/adefects-suite.sh frmtmb tests/testthat > dev/adefects-log/x.txt
#
# FRMTMB_BRMS_FIT_TESTS=true is the caller's to set for the gated tier.
set -e
PKG="$1"
DIR="$2"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
# NOT_CRAN is SET HERE, not left to the caller. Without it every
# skip_on_cran() block is skipped and a whole file reports
# pass=0 fail=0 err=0 skip=0, which reads exactly like a file with
# nothing in it. That happened to 43 core files on this lane's second
# pass, because the launcher did not export it.
export NOT_CRAN=true
n=0
for f in "$DIR"/test-*.R; do
  n=$((n + 1))
  "$R" dev/adefects-run.R "$PKG" "$f" 2>&1 |
    grep -E "^RESULT|^  BLOCK|^Error" || echo "RESULT $PKG $f NOOUTPUT"
done
echo "SUITE $PKG ran $n files"

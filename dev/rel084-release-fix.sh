#!/usr/bin/env bash
# The brief's merge steps 12 and 14 that are text edits: versions and
# floors in the DESCRIPTION files, and the release run scripts pointed at
# rellib-r23. Rerunnable: each edit reads "already applied" once done.
# CRLF is kept: every pattern stops before the end of its line. The NEWS
# sections are written by dev/rel084-news.R; codemeta.json is regenerated
# after this (dev/release/checklist.md, step 3).
#
#   bash dev/rel084-release-fix.sh
set -euo pipefail
REL=/c/Users/adf44/source/r/frmtmb-wt-rel084
cd "$REL"
ed1() {  # file, old, new
  if grep -qF -- "$3" "$1"; then echo "already applied: $1: $3"
  elif grep -qF -- "$2" "$1"; then
    perl -pi -e 'BEGIN { ($o, $n) = (shift, shift) } s/\Q$o\E/$n/' "$2" "$3" "$1"
    grep -qF -- "$3" "$1" || { echo "FAILED $1"; exit 1; }
    echo "applied: $1: $2 -> $3"
  else echo "NOT FOUND in $1: $2"; exit 1; fi
}
ed1 DESCRIPTION "Version: 0.83.0" "Version: 0.84.0"
# eam and survival fill the renamed cdf slot (lane cdfname): BREAKING, minor
ed1 extensions/frmtmb.eam/DESCRIPTION "Version: 0.13.2" "Version: 0.14.0"
ed1 extensions/frmtmb.eam/DESCRIPTION "frmtmb (>= 0.82.0)" "frmtmb (>= 0.84.0)"
ed1 extensions/frmtmb.survival/DESCRIPTION "Version: 0.3.0" "Version: 0.4.0"
ed1 extensions/frmtmb.survival/DESCRIPTION "frmtmb (>= 0.75.0)" "frmtmb (>= 0.84.0)"
# compat text only (a string in R/zzz.R): patch, floor unchanged
ed1 extensions/frmtmb.coupling/DESCRIPTION "Version: 0.6.1" "Version: 0.6.2"
ed1 dev/release/run-tests.R 'LIB <- "C:/Users/adf44/source/r/rellib-r22"' 'LIB <- "C:/Users/adf44/source/r/rellib-r23"'
ed1 dev/release/run-session.R 'LIB <- "C:/Users/adf44/source/r/rellib-r22"' 'LIB <- "C:/Users/adf44/source/r/rellib-r23"'
echo "RELEASE FIX DONE"

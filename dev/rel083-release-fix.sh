#!/usr/bin/env bash
# The brief's merge steps 12 and 14 that are text edits: versions and
# floors in the DESCRIPTION files, and the release run scripts pointed at
# rellib-r22. Rerunnable: each edit reads "already applied" once done.
# CRLF is kept: every pattern stops before the end of its line. The NEWS
# sections are written by dev/rel083-news.R; codemeta.json is regenerated
# after this (dev/release/checklist.md, step 3).
#
#   bash dev/rel083-release-fix.sh
set -euo pipefail
REL=/c/Users/adf44/source/r/frmtmb-wt-rel083
cd "$REL"
ed1() {  # file, old, new
  if grep -qF -- "$3" "$1"; then echo "already applied: $1: $3"
  elif grep -qF -- "$2" "$1"; then
    perl -pi -e 'BEGIN { ($o, $n) = (shift, shift) } s/\Q$o\E/$n/' "$2" "$3" "$1"
    grep -qF -- "$3" "$1" || { echo "FAILED $1"; exit 1; }
    echo "applied: $1: $2 -> $3"
  else echo "NOT FOUND in $1: $2"; exit 1; fi
}
ed1 DESCRIPTION "Version: 0.82.1" "Version: 0.83.0"
ed1 extensions/frmtmb.sample/DESCRIPTION "Version: 0.29.0" "Version: 0.30.0"
ed1 extensions/frmtmb.sample/DESCRIPTION "frmtmb (>= 0.82.0)" "frmtmb (>= 0.83.0)"
ed1 extensions/frmtmb.spatial/DESCRIPTION "Version: 0.4.0" "Version: 0.5.0"
ed1 extensions/frmtmb.spatial/DESCRIPTION "frmtmb (>= 0.82.0)" "frmtmb (>= 0.83.0)"
# tests only, but they need this core (quadseam2's and refits' tests)
ed1 extensions/frmtmb.latent/DESCRIPTION "Version: 0.7.4" "Version: 0.7.5"
ed1 extensions/frmtmb.latent/DESCRIPTION "frmtmb (>= 0.80.0)" "frmtmb (>= 0.83.0)"
ed1 dev/release/run-tests.R 'LIB <- "C:/Users/adf44/source/r/rellib-r21"' 'LIB <- "C:/Users/adf44/source/r/rellib-r22"'
ed1 dev/release/run-session.R 'LIB <- "C:/Users/adf44/source/r/rellib-r21"' 'LIB <- "C:/Users/adf44/source/r/rellib-r22"'
echo "RELEASE FIX DONE"

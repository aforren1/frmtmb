#!/usr/bin/env bash
# The brief's merge steps 10 and 12 that are text edits: versions and floors
# in the DESCRIPTION files, and the release run scripts pointed at
# rellib-r21. Rerunnable: each edit reads "already applied" once done.
# CRLF is kept: every pattern stops before the end of its line.
#
#   bash dev/rel082-release-fix.sh
set -euo pipefail
REL=/c/Users/adf44/source/r/frmtmb-wt-rel082
cd "$REL"
ed1() {  # file, old, new
  if grep -qF -- "$3" "$1"; then echo "already applied: $1: $3"
  elif grep -qF -- "$2" "$1"; then
    perl -pi -e 'BEGIN { ($o, $n) = (shift, shift) } s/\Q$o\E/$n/' "$2" "$3" "$1"
    grep -qF -- "$3" "$1" || { echo "FAILED $1"; exit 1; }
    echo "applied: $1: $2 -> $3"
  else echo "NOT FOUND in $1: $2"; exit 1; fi
}
ed1 DESCRIPTION "Version: 0.81.0" "Version: 0.82.0"
ed1 extensions/frmtmb.sample/DESCRIPTION "Version: 0.28.0" "Version: 0.29.0"
ed1 extensions/frmtmb.sample/DESCRIPTION "frmtmb (>= 0.81.0)" "frmtmb (>= 0.82.0)"
ed1 extensions/frmtmb.spatial/DESCRIPTION "Version: 0.3.0" "Version: 0.4.0"
ed1 extensions/frmtmb.spatial/DESCRIPTION "frmtmb (>= 0.76.0)" "frmtmb (>= 0.82.0)"
ed1 extensions/frmtmb.survival/DESCRIPTION "Version: 0.2.3" "Version: 0.3.0"
ed1 extensions/frmtmb.latent/DESCRIPTION "Version: 0.7.3" "Version: 0.7.4"
ed1 dev/release/run-tests.R 'LIB <- "C:/Users/adf44/source/r/rellib-r20"' 'LIB <- "C:/Users/adf44/source/r/rellib-r21"'
ed1 dev/release/run-session.R 'LIB <- "C:/Users/adf44/source/r/rellib-r20"' 'LIB <- "C:/Users/adf44/source/r/rellib-r21"'
# The release review's minor 4: frmtmb.eam and frmtmb.ode changed only
# tests this round, but those tests need 0.82.0, so each takes a patch
# with that floor and a NEWS section that says why.
ed1 extensions/frmtmb.eam/DESCRIPTION "Version: 0.13.1" "Version: 0.13.2"
ed1 extensions/frmtmb.eam/DESCRIPTION "frmtmb (>= 0.72.0)" "frmtmb (>= 0.82.0)"
ed1 extensions/frmtmb.ode/DESCRIPTION "Version: 0.7.3" "Version: 0.7.4"
ed1 extensions/frmtmb.ode/DESCRIPTION "frmtmb (>= 0.68.0)" "frmtmb (>= 0.82.0)"
news_top() {  # file, heading, text (one paragraph, already wrapped)
  if head -1 "$1" | tr -d '\r' | grep -qxF -- "$2"; then
    echo "already applied: $1: $2"
  else
    { printf '%s\r\n\r\n' "$2"; printf '%s\n' "$3" | sed 's/$/\r/'; printf '\r\n'; cat "$1"; } > "$1.new"
    mv "$1.new" "$1"; echo "applied: $1: $2"
  fi
}
news_top extensions/frmtmb.eam/NEWS.md "# frmtmb.eam 0.13.2" \
"Needs frmtmb 0.82.0: a test reads the message of frmtmb's standard-error
check for a parameter at the end of its range. No change to the
package's code."
news_top extensions/frmtmb.ode/NEWS.md "# frmtmb.ode 0.7.4" \
"Needs frmtmb 0.82.0: a test reads \`ranef()\` as one matrix per grouping
factor. No change to the package's code."
echo "RELEASE FIX DONE"

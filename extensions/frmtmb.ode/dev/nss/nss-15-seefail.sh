#!/bin/bash
# Every new test, run against a build with the thing it guards ABSENT.
#
# A test that pins a defect is worthless until it has been seen to fail,
# and "the symbol is missing" is not that. Each variant here keeps the
# whole API and removes one property of the correction, so the failure
# is behavioural:
#
#   none      the correction is not applied at all, which is bit for bit
#             what the base commit does (nss-09 section C)
#   joint     one ratio shared by every state, instead of one per state
#   nodamp    the noise floor removed from the denominator
#   hardcap   the stand-down replaced by a hard cap, C0 but not C1
#   outside   the correction built before the ss_extrapolate branch
#   oldreport the punch-round-0 diagnostic, built from the correction
#
# Script path: extensions/frmtmb.ode/dev/nss/nss-15-seefail.sh
set -u
SRC=C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode
WORK=C:/Users/adf44/source/r/nss-variants
LIBS=C:/Users/adf44/source/r/nsslib-variant
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="$LIBS;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"

rm -rf "$WORK" "$LIBS"
mkdir -p "$WORK" "$LIBS"

for V in none joint nodamp hardcap oldgate outside oldreport noundone; do
  rm -f "$SRC/dev/nss/nss-15-$V.log"
  rm -rf "$WORK/frmtmb.ode"
  cp -r "$SRC" "$WORK/frmtmb.ode"
  rm -rf "$WORK/frmtmb.ode/dev"
  if ! python "$SRC/dev/nss/nss-15-variant.py"         "$WORK/frmtmb.ode/R/ode.R" "$V"; then
    echo "=== variant: $V === PATCH FAILED, NOT RUN"
    continue
  fi
  # a variant that installed the shipped code is not a variant; the
  # first version of this harness did exactly that on all eight
  if cmp -s "$WORK/frmtmb.ode/R/ode.R" "$SRC/R/ode.R"; then
    echo "=== variant: $V === SOURCE UNCHANGED, NOT RUN"
    continue
  fi
  "$R/R.exe" CMD INSTALL --library="$LIBS" "$WORK/frmtmb.ode" \
      > "$SRC/dev/nss/nss-15-install-$V.log" 2>&1
  echo "=== variant: $V ==="
  R_LIBS="$LIBS;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6" \
    "$R/Rscript.exe" "$SRC/dev/nss/nss-15-run.R" \
    > "$SRC/dev/nss/nss-15-$V.log" 2>&1
  grep -E "^(FILE|── [0-9]+\.)" "$SRC/dev/nss/nss-15-$V.log"
done

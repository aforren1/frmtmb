#!/bin/bash
# One see-it-fail variant, with its own work directory and library so
# that two runs cannot tread on each other. The whole-set driver
# `nss-15-seefail.sh` deadlocked twice when it was started more than
# once by accident, because both copies shared `$WORK` and `$LIBS`.
#
# Usage: nss-15-one.sh <variant>
# Script path: extensions/frmtmb.ode/dev/nss/nss-15-one.sh
set -u
V=${1:?variant}
SRC=C:/Users/adf44/source/r/frmtmb-wt-nss/extensions/frmtmb.ode
WORK=C:/Users/adf44/source/r/nss-var-$V
LIBS=C:/Users/adf44/source/r/nsslib-var-$V
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="$LIBS;C:/Users/adf44/source/r/rellib-0552;C:/Users/adf44/source/r/pinlib;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true

rm -rf "$WORK" "$LIBS" "$SRC/dev/nss/nss-15-$V.log"
mkdir -p "$WORK" "$LIBS"
cp -r "$SRC" "$WORK/frmtmb.ode"
rm -rf "$WORK/frmtmb.ode/dev"
if ! python "$SRC/dev/nss/nss-15-variant.py" \
      "$WORK/frmtmb.ode/R/ode.R" "$V"; then
  echo "$V: PATCH FAILED, NOT RUN"
  exit 1
fi
if cmp -s "$WORK/frmtmb.ode/R/ode.R" "$SRC/R/ode.R"; then
  echo "$V: SOURCE UNCHANGED, NOT RUN"
  exit 1
fi
"$R/R.exe" CMD INSTALL --library="$LIBS" "$WORK/frmtmb.ode" \
    > "$SRC/dev/nss/nss-15-install-$V.log" 2>&1 || {
  echo "$V: INSTALL FAILED"; exit 1; }
"$R/Rscript.exe" "$SRC/dev/nss/nss-15-run.R" \
    > "$SRC/dev/nss/nss-15-$V.log" 2>&1
grep -E "^FILE|^BROKE" "$SRC/dev/nss/nss-15-$V.log" ||
  echo "$V: NO SUMMARY LINE"
rm -rf "$WORK" "$LIBS"

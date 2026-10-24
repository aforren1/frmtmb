#!/usr/bin/env bash
# Does frmtmb.survival need this round's core? It is installed ALONE into
# a scratch library over rellib-r11 (frmtmb 0.73.0), from a copy whose
# DESCRIPTION floor is lowered so that the install does not refuse on
# the version, then loaded and asked for a cox() fit and a curve off a
# royston_parmar() fit. Adapted from dev/rel073-floor.sh.
#
#   bash dev/splitsurv-floor.sh <scratch dir outside the worktree>
set -u
S="$1"
ROOT=C:/Users/adf44/source/r/frmtmb-wt-splitsurv
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
FL="$S/floor-lib"; SRC="$S/floor-src/frmtmb.survival"
rm -rf "$FL" "$S/floor-src"; mkdir -p "$FL" "$S/floor-src"
cp -r "$ROOT/extensions/frmtmb.survival" "$SRC"
sed -i -e 's/frmtmb (>= 0.73.1)/frmtmb (>= 0.61.0)/' "$SRC/DESCRIPTION"
grep -n 'frmtmb (>=' "$SRC/DESCRIPTION"
echo "== install over rellib-r11"
R_LIBS="$FL;C:/Users/adf44/source/r/rellib-r11;C:/Users/adf44/AppData/Local/R/win-library/4.6" \
  "$R" CMD INSTALL --library="$FL" --no-multiarch "$SRC" 2>&1 |
  grep -E 'DONE|ERROR|Error|error|not exported|installation of package'
echo "== load and fit"
"$RS" "$ROOT/dev/splitsurv-floor-probe.R" "$FL" 2>&1 | tail -20

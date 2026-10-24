#!/usr/bin/env bash
# Roxygenise and install the lane's packages into its private library,
# core first. Usage: bash dev/splitsurv-install.sh [pkg ...]
# where pkg is frmtmb or an extension name; no argument installs core,
# frmtmb.survival, frmtmb.sample and frmtmb.learn.
set -eu
ROOT=C:/Users/adf44/source/r/frmtmb-wt-splitsurv
LIB=C:/Users/adf44/source/r/wt-splitsurv-lib
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r11;C:/Users/adf44/AppData/Local/R/win-library/4.6"
cd "$ROOT"
pk=("$@")
[ ${#pk[@]} -gt 0 ] || pk=(frmtmb frmtmb.survival frmtmb.sample frmtmb.learn)
for p in "${pk[@]}"; do
  if [ "$p" = frmtmb ]; then d="$ROOT"; else d="$ROOT/extensions/$p"; fi
  echo "== roxygenise $p"
  "$RS" -e "roxygen2::roxygenise('$d')" 2>&1 | grep -v "^ℹ\|Loading\|Writing '" | tail -20
  echo "== install $p"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$d" 2>&1 | tail -3
done

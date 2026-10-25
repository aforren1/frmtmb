#!/usr/bin/env bash
# Review: spot-check two floors of dev/rel075-floor.sh (frmtmb.survival,
# frmtmb.spatial) and frmtmb.latent's test-mixture-mvn.R, with every
# scratch library under rev-rel075-lib and Rtools NOT first on PATH.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel075
LIB=C:/Users/adf44/source/r/rev-rel075-lib
R12=C:/Users/adf44/source/r/rellib-r12
R13=C:/Users/adf44/source/r/rellib-r13
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mkdir -p /c/Users/adf44/source/r/rev-rel075-lib/src
cd "$ROOT" || exit 1
load_line() {
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
inst_line() {
  R_LIBS="$1;$2;$USR" "$R" CMD INSTALL --library="$1" --no-multiarch "$3" 2>&1 |
    grep -E 'DONE|ERROR|Error|error|not exported|installation of package|unused argument' | sed 's/^/  /' | head -8
}
for p in frmtmb.survival frmtmb.spatial frmtmb.latent; do
  echo "== $p $(grep '^Version:' extensions/$p/DESCRIPTION) floor $(grep -o 'frmtmb (>= [0-9.]*)' extensions/$p/DESCRIPTION)"
  SRC=$LIB/src/$p; rm -rf "$SRC"; cp -r "extensions/$p" "$SRC"
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$SRC/DESCRIPTION"
  FL=$LIB/old-$p; rm -rf "$FL"; mkdir -p "$FL"
  echo "-- old (rellib-r12, floor lowered)"; inst_line "$FL" "$R12" "$SRC"
  [ -d "$FL/$p" ] && load_line "$FL" "$R12" "$p"
  FN=$LIB/new-$p; rm -rf "$FN"; mkdir -p "$FN"
  echo "-- new (rellib-r13)"; inst_line "$FN" "$R13" "extensions/$p"
  load_line "$FN" "$R13" "$p"
done
echo "-- frmtmb.latent test-mixture-mvn.R on old core"
NOT_CRAN=true ROUND_BASE=$R12 "$RS" dev/ciharden-run1.R "$LIB/old-frmtmb.latent" frmtmb.latent \
  extensions/frmtmb.latent/tests/testthat/test-mixture-mvn.R 2>&1 | grep -E '^lib:|^RESULT|^DETAIL.*(FAIL|ERROR)' | cut -c1-300
echo "-- frmtmb.latent test-mixture-mvn.R on new core"
NOT_CRAN=true ROUND_BASE=$R13 "$RS" dev/ciharden-run1.R "$LIB/new-frmtmb.latent" frmtmb.latent \
  extensions/frmtmb.latent/tests/testthat/test-mixture-mvn.R 2>&1 | grep -E '^lib:|^RESULT|^DETAIL.*(FAIL|ERROR)' | cut -c1-300
rm -rf "$LIB/src"
echo FLOOR-REVIEW DONE

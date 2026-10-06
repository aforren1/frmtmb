#!/usr/bin/env bash
# Install the release worktree into rellib-r6: roxygenise core, install
# it, then roxygenise every extension against it and install them in
# dependency order (frmtmb.learn imports frmtmb.eam). drmTMB 0.7.0 is
# copied once from rellib-r5, which holds the CRAN build, so the gated
# drmTMB file does not skip.
#
#   bash dev/rel068-install.sh            # core and all seven extensions
#   bash dev/rel068-install.sh core       # core only, no roxygen
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
REL=C:/Users/adf44/source/r/rellib-r6
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$REL;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
mkdir -p "$REL"
if [ ! -d "$REL/drmTMB" ]; then
  cp -r C:/Users/adf44/source/r/rellib-r5/drmTMB "$REL/"
fi
inst() {
  echo "== INSTALL $1 $(date +%T)"
  "$R" CMD INSTALL --library="$REL" --no-multiarch "$2" 2>&1 |
    grep -E 'DONE|ERROR|Error|WARN' | head -8
}
if [ "${1:-all}" = core ]; then
  inst frmtmb "$ROOT"
  exit 0
fi
"$RS" "$ROOT/dev/rel068-rox.R" frmtmb 2>&1 | grep -vE '^(Writing|Loading|ℹ|i Loading)' | tail -15
inst frmtmb "$ROOT"
"$RS" "$ROOT/dev/rel068-rox.R" frmtmb.eam frmtmb.sample frmtmb.spline frmtmb.latent \
  frmtmb.learn frmtmb.coupling frmtmb.ode 2>&1 | grep -E '^== rox|ROX DONE|Error|Warn|Writing'
for p in frmtmb.eam frmtmb.latent frmtmb.ode frmtmb.spline frmtmb.coupling \
         frmtmb.sample frmtmb.learn; do
  inst "$p" "$ROOT/extensions/$p"
done
"$RS" -e '.libPaths(c("C:/Users/adf44/source/r/rellib-r6", "C:/Users/adf44/AppData/Local/R/win-library/4.6")); for (p in c("frmtmb", paste0("frmtmb.", c("eam","latent","ode","spline","coupling","sample","learn")), "drmTMB")) cat(sprintf("%-16s %s %s\n", p, as.character(packageVersion(p)), dirname(find.package(p))))'
echo "INSTALL-RELEASE DONE"

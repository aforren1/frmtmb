#!/usr/bin/env bash
# Install the release worktree into a library: roxygenise core, install
# it, then roxygenise every extension against it and install them in
# dependency order. frmtmb.learn imports frmtmb.eam; frmtmb.sample
# suggests eam, latent, ode and survival, so they go first. frmtmb.spline
# is retired and frmtmb.survival and frmtmb.spatial are new. drmTMB is
# copied once from rellib-r11, so the gated drmTMB file does not skip.
# Adapted from dev/rel073-install.sh.
#
#   bash dev/rel074-install.sh <lib>          # core and all eight
#   bash dev/rel074-install.sh <lib> core     # core only, no roxygen
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
LIB=${1:?library}
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export REL074_LIB="$LIB"
case "$LIB" in *rellib-r11*|*rellib-r10*|*rellib-r9*)
  echo "rellib-r9, r10 and r11 are read-only"; exit 1;; esac
mkdir -p "$LIB"
if [ ! -d "$LIB/drmTMB" ]; then
  cp -r C:/Users/adf44/source/r/rellib-r11/drmTMB "$LIB/"
fi
inst() {
  echo "== INSTALL $1 $(date +%T)"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$2" 2>&1 |
    grep -E 'DONE|ERROR|Error|WARN' | head -8
}
EXT="frmtmb.eam frmtmb.latent frmtmb.ode frmtmb.survival frmtmb.spatial
     frmtmb.coupling frmtmb.sample frmtmb.learn"
if [ "${2:-all}" = core ]; then
  inst frmtmb "$ROOT"
  exit 0
fi
"$RS" "$ROOT/dev/rel074-rox.R" frmtmb 2>&1 | grep -vE '^(Writing|Loading|ℹ|i Loading)' | tail -15
inst frmtmb "$ROOT"
"$RS" "$ROOT/dev/rel074-rox.R" $EXT 2>&1 | grep -E '^== rox|ROX DONE|Error|Warn|Writing|Deleting'
for p in $EXT; do
  inst "$p" "$ROOT/extensions/$p"
done
"$RS" -e '.libPaths(c(Sys.getenv("REL074_LIB"), "C:/Users/adf44/AppData/Local/R/win-library/4.6")); for (p in c("frmtmb", paste0("frmtmb.", c("eam","latent","ode","survival","spatial","coupling","sample","learn")), "drmTMB")) cat(sprintf("%-16s %s %s\n", p, as.character(packageVersion(p)), dirname(find.package(p))))'
echo "INSTALL DONE"

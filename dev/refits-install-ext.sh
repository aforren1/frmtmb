#!/usr/bin/env bash
# Install every extension from the worktree into the lane library, after
# core, so each loads the lane's core. No roxygenise: their sources are
# unchanged.
WT=C:/Users/adf44/source/r/frmtmb-wt-refits
LIB=C:/Users/adf44/source/r/wt-refits-lib
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
for e in latent survival spatial ode coupling sample eam learn; do
  echo "== frmtmb.$e"
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$WT/extensions/frmtmb.$e" 2>&1 | tail -2
done

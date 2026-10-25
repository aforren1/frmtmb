#!/usr/bin/env bash
# Roxygenise and install core and the named extensions from this
# worktree into the lane library, in that order.
#   bash dev/smallfix-install.sh [frmtmb.latent frmtmb.survival ...]
# Log: dev/smallfix-out/install.txt
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-smallfix
LIB=C:/Users/adf44/source/r/wt-smallfix-lib
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
log=$WT/dev/smallfix-out/install.txt
: > "$log"
one() {
  local dir=$1
  echo "== $dir" >> "$log"
  "$R/Rscript.exe" -e "roxygen2::roxygenise('$dir')" >> "$log" 2>&1 || return 1
  "$R/R.exe" CMD INSTALL --library="$LIB" --no-multiarch "$dir" >> "$log" 2>&1 || return 1
}
one "$WT" || { echo "core FAILED"; exit 1; }
for p in "$@"; do
  one "$WT/extensions/$p" || { echo "$p FAILED"; exit 1; }
done
grep -E "^\* DONE|^== |ERROR|Warning" "$log"

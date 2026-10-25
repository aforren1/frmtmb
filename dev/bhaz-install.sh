#!/usr/bin/env bash
# Roxygenise and install core and frmtmb.survival from this worktree into
# the bhaz lane library. Usage: bash dev/bhaz-install.sh [core] [survival]
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-bhaz
LIB=C:/Users/adf44/source/r/wt-bhaz-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_LIBS="$LIB;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
mkdir -p "$WT/dev/bhaz-log"
for what in "$@"; do
  case "$what" in
    core) dir="$WT" ;;
    survival) dir="$WT/extensions/frmtmb.survival" ;;
    *) dir="$WT/extensions/$what" ;;
  esac
  "$RS" -e ".libPaths(strsplit(Sys.getenv('R_LIBS'), ';')[[1]]); roxygen2::roxygenise('$dir')" \
    > "$WT/dev/bhaz-log/rox-$what.txt" 2>&1
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$dir" \
    > "$WT/dev/bhaz-log/install-$what.txt" 2>&1
  tail -2 "$WT/dev/bhaz-log/install-$what.txt"
done

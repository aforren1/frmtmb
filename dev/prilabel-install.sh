#!/usr/bin/env bash
# Roxygenise and install lane prilabel's core (and any extension named
# as an argument, e.g. frmtmb.sample) into the lane's private library.
set -eu
WT=C:/Users/adf44/source/r/frmtmb-wt-prilabel
LIB=C:/Users/adf44/source/r/wt-prilabel-lib
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
R9=C:/Users/adf44/source/r/rellib-r9
UL=C:/Users/adf44/AppData/Local/R/win-library/4.6
export R_LIBS="$LIB;$R9;$UL"
if [ "$#" -eq 0 ]; then set -- core; fi
for p in "$@"; do
  if [ "$p" = core ]; then d="$WT"; else d="$WT/extensions/$p"; fi
  "$RS" -e "roxygen2::roxygenise('$d')" > "$WT/dev/prilabel-log/rox-$p.txt" 2>&1
  "$R" CMD INSTALL --library="$LIB" --no-multiarch "$d" \
    > "$WT/dev/prilabel-log/install-$p.txt" 2>&1
  tail -2 "$WT/dev/prilabel-log/install-$p.txt"
done

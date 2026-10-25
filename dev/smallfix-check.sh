#!/usr/bin/env bash
# R CMD build and R CMD check --as-cran of one package of this
# worktree, inside dev/smallfix-check/<package>/, with NOT_CRAN=true,
# pandoc and TinyTeX on PATH, and R_LIBS = the lane library, then
# rellib-r12, then the user library (lane-env.md, Deliverables 5).
#   bash dev/smallfix-check.sh <. | extensions/frmtmb.xxx>
# Log: dev/smallfix-check/<package>/check.txt; Status line at the end.
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-smallfix
src=$WT/$1
pkg=$(sed -n 's/^Package: //p' "$src/DESCRIPTION" | tr -d '\r')
dir=$WT/dev/smallfix-check/$pkg
rm -rf "$dir"; mkdir -p "$dir"
cd "$dir" || exit 1
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-smallfix-lib;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
"$R" CMD build "$src" > build.txt 2>&1
tgz=$(ls "$pkg"_*.tar.gz 2>/dev/null | head -1)
[ -n "$tgz" ] || { echo "$pkg: build failed"; tail -5 build.txt; exit 1; }
"$R" CMD check --as-cran "$tgz" > check.txt 2>&1
echo "$pkg: $(grep '^Status' check.txt)"

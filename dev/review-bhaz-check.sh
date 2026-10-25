#!/usr/bin/env bash
# Reviewer: R CMD build + check --as-cran of frmtmb.survival in
# dev/review-bhaz-check/.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-bhaz
src="$WT/extensions/frmtmb.survival"
out="$WT/dev/review-bhaz-check"; rm -rf "$out"; mkdir -p "$out"; cd "$out"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-bhaz-lib;C:/Users/adf44/source/r/rellib-r12;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
"$R" CMD build "$(cygpath -w "$src")" > build.txt 2>&1
tb=$(ls *.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > check.txt 2>&1
grep "^Status" check.txt

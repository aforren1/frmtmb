#!/bin/bash
# R CMD build and check --as-cran of one package of the lane, inside
# dev/reslot-check/<name>. Usage: bash dev/reslot-check.sh <dir>
# where <dir> is "." for core or extensions/<pkg>.
WT=C:/Users/adf44/source/r/frmtmb-wt-reslot
src="$WT/$1"; nm=$(basename "$(cd "$src" && pwd)")
[ "$1" = "." ] && nm=frmtmb
OUT="$WT/dev/reslot-check/$nm"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
export R_LIBS="C:/Users/adf44/source/r/wt-reslot-lib;C:/Users/adf44/source/r/rellib-r14;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$src" > build.log 2>&1
tgz=$(ls -t *.tar.gz 2>/dev/null | head -1)
[ -n "$tgz" ] || { echo "NO TARBALL"; exit 1; }
"$R" CMD check --as-cran "$tgz" > check.log 2>&1
grep -a "^Status" check.log

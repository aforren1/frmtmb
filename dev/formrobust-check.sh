#!/usr/bin/env bash
# R CMD build + check --as-cran of one package of the worktree, inside
# dev/formrobust-check/<name>. Usage: bash dev/formrobust-check.sh <path>
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-formrobust
src="$1"; nm=$(basename "$src")
[ "$nm" = "frmtmb-wt-formrobust" ] && nm=frmtmb
OUT="$WT/dev/formrobust-check/$nm"
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="C:/Users/adf44/source/r/wt-formrobust-lib;C:/Users/adf44/source/r/rellib-r4;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
"$R" CMD build "$src" > build.log 2>&1
tgz=$(ls -t *.tar.gz | head -1)
"$R" CMD check --as-cran "$tgz" > check.log 2>&1
grep -a "^Status" check.log

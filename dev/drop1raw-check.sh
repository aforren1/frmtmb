#!/bin/sh
# R CMD build and check --as-cran of one package (core: "."), inside
# dev/drop1raw-check/<name>. Log: dev/drop1raw-check/<name>/check.txt
src="$1"; name="$2"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-drop1raw
D="$ROOT/dev/drop1raw-check/$name"
rm -rf "$D"; mkdir -p "$D"; cd "$D" || exit 1
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_LIBS="C:/Users/adf44/source/r/wt-drop1raw-lib;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$ROOT/$src" > build.txt 2>&1 || { tail build.txt; exit 1; }
tgz=$(ls *.tar.gz | head -1)
"$R" CMD check --as-cran "$tgz" > check.txt 2>&1
grep "^Status" check.txt

#!/bin/bash
# Lane meme: R CMD check --as-cran of one package, built and checked in
# dev/meme-check/. Usage: bash dev/meme-check.sh <source dir> <name>
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-meme
SRC="$1"; NM="$2"
OUT=$WT/dev/meme-check/$NM
rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
export R_LIBS="C:/Users/adf44/source/r/wt-meme-lib;C:/Users/adf44/source/r/rellib-r15;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-meme/dev/stan-cache
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$SRC" > build.txt 2>&1
TB=$(ls *.tar.gz | head -1)
"$R" CMD check --as-cran "$TB" > check.txt 2>&1
grep "^Status" check.txt

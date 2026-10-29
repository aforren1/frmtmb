#!/bin/sh
# Lane rescorpp: R CMD build and check --as-cran of one package, inside
# dev/rescorpp-check/<pkg>. Usage: sh dev/rescorpp-check.sh core|sample
WT=C:/Users/adf44/source/r/frmtmb-wt-rescorpp
R="/c/Program Files/R/R-4.6.1/bin"
case $1 in core) src=$WT ;; sample) src=$WT/extensions/frmtmb.sample ;; esac
out=$WT/dev/rescorpp-check/$1
rm -rf "$out"; mkdir -p "$out"; cd "$out" || exit 1
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-rescorpp-lib;C:/Users/adf44/source/r/rellib-r15;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
"$R/R.exe" CMD build "$src" > build.txt 2>&1 || { tail build.txt; exit 1; }
tb=$(ls *.tar.gz)
"$R/R.exe" CMD check --as-cran "$tb" > check.txt 2>&1
grep -n "^Status" -A 3 check.txt

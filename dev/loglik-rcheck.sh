#!/bin/sh
# Lane loglik: R CMD build and R CMD check --as-cran of one package of
# this tree, inside dev/loglik-check/<pkg>/.
#   sh dev/loglik-rcheck.sh <dir relative to the tree> <package name>
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-loglik
R="/c/Program Files/R/R-4.6.1/bin"
export R_LIBS="C:/Users/adf44/source/r/wt-loglik-lib;C:/Users/adf44/source/r/rellib-r13;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_LIBS_USER="$R_LIBS"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
out="$WT/dev/loglik-check/$2"
rm -rf "$out"; mkdir -p "$out"; cd "$out"
"$R/R.exe" CMD build "$WT/$1"
tb=$(ls "$2"_*.tar.gz)
"$R/R.exe" CMD check --as-cran "$tb"

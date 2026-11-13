#!/usr/bin/env bash
# Lane sampfront: R CMD build and R CMD check --as-cran of one package,
# inside dev/sampfront-check/<package>/, against the lane library, then
# rellib-r20, then the user library.
#   bash dev/sampfront-check.sh [frmtmb.sample|frmtmb] \
#     > dev/sampfront-out/check-<package>.txt 2>&1
set -u
P="${1:-frmtmb.sample}"
ROOT="C:/Users/adf44/source/r/frmtmb-wt-sampfront"
SRC="$ROOT/extensions/$P"
[ "$P" = frmtmb ] && SRC="$ROOT"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
PD="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TT=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
MG=$RT/x86_64-w64-mingw32.static.posix/bin
export PATH="$PD:$TT:$PATH:$RT/usr/bin:$MG"
export RSTUDIO_PANDOC="C:${PD#/c}"
RR=C:/Users/adf44/source/r
UL=C:/Users/adf44/AppData/Local/R/win-library/4.6
export R_LIBS="$RR/wt-sampfront-lib;$RR/rellib-r20;$UL"
export R_LIBS_USER=""
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
CK="$ROOT/dev/sampfront-check/$P"
rm -rf "$CK"; mkdir -p "$CK"; cd "$CK" || exit 1
"$R" CMD build "$SRC"
tb=$(ls "${P}"_*.tar.gz)
"$R" CMD check --as-cran "$tb"
grep "^Status" "$P.Rcheck/00check.log"

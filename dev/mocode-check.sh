#!/usr/bin/env bash
# R CMD build and check --as-cran of one package, inside
# dev/mocode-check/<pkg>/, against the lane library.
#   bash dev/mocode-check.sh frmtmb|frmtmb.sample
set -eu
pkg="$1"
W=/c/Users/adf44/source/r/frmtmb-wt-mocode
if [ "$pkg" = frmtmb ]; then src="$W"; else src="$W/extensions/$pkg"; fi
D="$W/dev/mocode-check/$pkg"
rm -rf "$D"; mkdir -p "$D"; cd "$D"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RT=/c/rtools45
PD="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TT=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PD:$TT:$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-mocode-lib;C:/Users/adf44/source/r/rellib-r13;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$W/dev/stan-cache"
unset FRMTMB_BRMS_FIT_TESTS
"$R" CMD build "$src" > build.txt 2>&1
tgz=$(ls "$pkg"_*.tar.gz)
"$R" CMD check --as-cran "$tgz" > check.txt 2>&1 || true
grep -E "^Status|WARNING|ERROR|NOTE" check.txt

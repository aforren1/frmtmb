#!/usr/bin/env bash
# Lane ndaterm: R CMD build and check --as-cran of one package, inside
# dev/ndaterm-check/<pkg>/. Usage: bash dev/ndaterm-check.sh <pkg>
set -eu
p="$1"
WT=/c/Users/adf44/source/r/frmtmb-wt-ndaterm
src="$WT"; [ "$p" = frmtmb ] || src="$WT/extensions/$p"
out="$WT/dev/ndaterm-check/$p"
rm -rf "$out"; mkdir -p "$out"; cd "$out"
R="/c/Program Files/R/R-4.6.1/bin/R"
export R_LIBS="C:/Users/adf44/source/r/wt-ndaterm-lib;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
"$R" CMD build "$src" > build.log 2>&1
tgz=$(ls "$p"_*.tar.gz)
"$R" CMD check --as-cran "$tgz" > check.log 2>&1 || true
grep -E "^Status|WARNING|ERROR|NOTE" check.log | tail -15

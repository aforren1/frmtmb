#!/usr/bin/env bash
# Lane priors2: R CMD build and R CMD check --as-cran of one package of
# the worktree, inside dev/priors2-check/<pkg>, with the lane library,
# then rellib-r21, then the user library.
#   bash dev/priors2-check.sh frmtmb | frmtmb.sample | frmtmb.spatial
set -eu
p="$1"
WT=C:/Users/adf44/source/r/frmtmb-wt-priors2
src="$WT"; [ "$p" = frmtmb ] || src="$WT/extensions/$p"
R="/c/Program Files/R/R-4.6.1/bin"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-priors2-lib;C:/Users/adf44/source/r/rellib-r21;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
out="$WT/dev/priors2-check/$p"
rm -rf "$out"; mkdir -p "$out"; cd "$out"
"$R/R.exe" CMD build "$src" > build.log 2>&1 || { echo "build failed"; tail -5 build.log; exit 1; }
tb=$(ls "$p"_*.tar.gz)
"$R/R.exe" CMD check --as-cran "$tb" > check.log 2>&1 || true
grep -aE "^Status|WARNING|ERROR|NOTE" check.log | head -20

#!/usr/bin/env bash
# Lane adhess: R CMD build and R CMD check --as-cran of one package of
# the worktree, inside dev/adhess-check/<pkg>/.
#   bash dev/adhess-check.sh <path relative to the worktree, "." for core>
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-adhess
src="$W/$1"
name=$(basename "$(cd "$src" && pwd)")
[ "$1" = "." ] && name=frmtmb
D="$W/dev/adhess-check/$name"
rm -rf "$D"; mkdir -p "$D"; cd "$D"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PANDOC:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export R_LIBS="C:/Users/adf44/source/r/wt-adhess-lib;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$W/dev/stan-cache"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$src" > build.txt 2>&1
tb=$(ls ${name}_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > check.txt 2>&1
grep -E "^Status|WARNING|ERROR|NOTE" check.txt

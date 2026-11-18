#!/usr/bin/env bash
# Lane gpselimit: R CMD build and R CMD check --as-cran of frmtmb.sample
# from this worktree (one test file changed), inside
# dev/gpselimit-check/frmtmb.sample/, with the lane's library first.
# Usage: bash dev/gpselimit-check-sample.sh
set -u
WT=C:/Users/adf44/source/r/frmtmb-wt-gpselimit
CK=/c/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-check/frmtmb.sample
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
PANDOC="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
TEX=/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows
export PATH="$PANDOC:$TEX:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export R_LIBS="C:/Users/adf44/source/r/wt-gpselimit-lib;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export NOT_CRAN=true
export FRMTMB_STAN_CACHE=$WT/dev/stan-cache
rm -rf "$CK"; mkdir -p "$CK"; cd "$CK"
"$R" CMD build "$WT/extensions/frmtmb.sample" > build.txt 2>&1 ||
  { tail build.txt; exit 1; }
tb=$(ls frmtmb.sample_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tb" > check.txt 2>&1
grep -a "^Status:" check.txt
grep -a -A3 "WARNING\|ERROR\|NOTE" check.txt | head -60

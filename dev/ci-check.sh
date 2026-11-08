#!/usr/bin/env bash
# R CMD check --as-cran on core from this worktree, once, at the end of
# lane ci, with the release driver's settings (dev/release/run-check.ps1)
# and the lane library first.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-ci
OUT=$WT/dev/ci-check
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="C:/Users/adf44/source/r/wt-ci-lib;C:/Users/adf44/source/r/rellib-r18;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export R_PARALLEL_PORT=11841
rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT"
"$R" CMD build "$WT" > build.log 2>&1
tgz=$(ls -t frmtmb_*.tar.gz 2>/dev/null | head -1)
[ -n "$tgz" ] || { echo "NO TARBALL BUILT"; exit 1; }
"$R" CMD check --as-cran "$tgz" > check.log 2>&1
grep -E "^Status:" check.log

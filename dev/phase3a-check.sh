#!/usr/bin/env bash
# R CMD build (with vignettes) and R CMD check --as-cran for one
# extension, in a directory outside the worktree. The check's
# testthat.Rout is copied into dev/phase3a-log/ before anything else,
# because a passing check deletes nothing but a later cleanup would.
# Usage: phase3a-check.sh <pkg>
set -u
pkg="$1"
W=C:/Users/adf44/source/r/frmtmb-wt-phase3a
C=C:/Users/adf44/source/r/phase3a-check/$pkg
LOG=$W/dev/phase3a-log
rm -rf "$C"; mkdir -p "$C"; cd "$C"
export R_LIBS="C:/Users/adf44/source/r/phase3a-lib;C:/Users/adf44/source/r/rellib-r3;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export R_LIBS_USER="C:/Users/adf44/source/r/phase3a-lib"
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH"
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export NOT_CRAN=true
R="C:/Program Files/R/R-4.6.1/bin/R.exe"
"$R" CMD build "$W/extensions/$pkg" > "$LOG/check-$pkg-build.log" 2>&1
tgz=$(ls "$pkg"_*.tar.gz | head -1)
"$R" CMD check --as-cran "$tgz" > "$LOG/check-$pkg.log" 2>&1
cp "$pkg.Rcheck/00check.log" "$LOG/check-$pkg-00check.log" 2>/dev/null
cp "$pkg.Rcheck/tests/testthat.Rout" "$LOG/check-$pkg-testthat.Rout" 2>/dev/null
cp "$pkg.Rcheck/tests/testthat.Rout.fail" "$LOG/check-$pkg-testthat.Rout.fail" 2>/dev/null
grep "^Status" "$LOG/check-$pkg.log"

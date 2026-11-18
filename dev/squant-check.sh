#!/usr/bin/env bash
# Lane squant: R CMD build and R CMD check --as-cran of core, adapted
# from dev/rel083-check.sh. Built with vignettes and checked with the
# manual; scratch in dev/squant-check/frmtmb/, the check log is the
# evidence. The tarball is passed by its absolute path.
# Run: bash dev/squant-check.sh
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-squant
OUT=$ROOT/dev/squant-check
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="C:/Users/adf44/source/r/wt-squant-lib;C:/Users/adf44/source/r/rellib-r22;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
: "${TMP:?TMP is unset}"
stray=$(cd "$ROOT" && find . -path ./dev -prune -o \( -name "testthat-problems.rds" \
  -o -name "_problems" -o -name "Rplots.pdf" -o -name "*.Rcheck" \) -print)
if [ -n "$stray" ]; then echo "STRAY FILES, not building:"; echo "$stray"; exit 1; fi
d=$OUT/frmtmb; rm -rf "$d"; mkdir -p "$d"; cd "$d" || exit 1
echo "== frmtmb build $(date +%T)"
"$R" CMD build "$ROOT" > build.log 2>&1
tgz=$(ls -t frmtmb_*.tar.gz 2>/dev/null | head -1)
if [ -z "$tgz" ]; then echo "NO TARBALL"; tail -5 build.log; exit 1; fi
echo "== frmtmb check $tgz $(date +%T)"
"$R" CMD check --as-cran "$d/$tgz" > check.log 2>&1
grep -E "^Status:" check.log || { echo "NO STATUS"; tail -5 check.log; }
echo "CHECK DONE $(date +%T)"

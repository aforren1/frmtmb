#!/usr/bin/env bash
# R CMD build and R CMD check --as-cran, one package at a time, for the
# packages whose version changed at 0.69.0. Built WITH vignettes and
# checked WITH the manual (dev/release/run-check.ps1 says why). Scratch
# lives in dev/release-check-069/<pkg>/; the check logs are the evidence.
set -u
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
OUT=$ROOT/dev/release-check-069
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
export R_LIBS="C:/Users/adf44/source/r/rellib-r7;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:\\Program Files\\RStudio\\resources\\app\\bin\\quarto\\bin\\tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
: "${TMP:?TMP is unset}"
# a tree with test leftovers is not clean: R CMD build would ship them
# (the surface review's process note: testthat-problems.rds is not in
# .Rbuildignore)
stray=$(cd "$ROOT" && find . -path ./dev -prune -o \( -name "testthat-problems.rds" \
  -o -name "_problems" -o -name "Rplots.pdf" -o -name "*.Rcheck" \) -print)
if [ -n "$stray" ]; then echo "STRAY FILES, not building:"; echo "$stray"; exit 1; fi
for p in "$@"; do
  src=$ROOT; [ "$p" = frmtmb ] || src=$ROOT/extensions/$p
  d=$OUT/$p; rm -rf "$d"; mkdir -p "$d"; cd "$d" || exit 1
  echo "== $p build $(date +%T)"
  "$R" CMD build "$src" > build.log 2>&1
  tgz=$(ls -t ${p}_*.tar.gz 2>/dev/null | head -1)
  if [ -z "$tgz" ]; then echo "NO TARBALL for $p"; tail -5 build.log; continue; fi
  echo "== $p check $tgz $(date +%T)"
  "$R" CMD check --as-cran "$tgz" > check.log 2>&1
  grep -E "^Status:" check.log || { echo "NO STATUS for $p"; tail -5 check.log; }
done
echo "CHECK DONE $(date +%T)"

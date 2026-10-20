#!/usr/bin/env bash
# R CMD check --as-cran on the packages lane testiso changed, built with
# vignettes inside dev/testiso-check/ (dev/release/run-check.ps1 settings)
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
OUT="$ROOT/dev/testiso-check"
export R_LIBS="C:/Users/adf44/source/r/wt-testiso-lib;C:/Users/adf44/source/r/rellib-r9;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
: "${TMP:?TMP is unset}"
mkdir -p "$OUT"; cd "$OUT"
# the packages named (core as "frmtmb"); by default the three of round 0
[ $# -gt 0 ] || set -- frmtmb frmtmb.ode frmtmb.sample
for n in "$@"; do
  if [ "$n" = frmtmb ]; then p="$ROOT"; else p="$ROOT/extensions/$n"; fi
  rm -rf "$n.Rcheck" ${n}_*.tar.gz
  "$R" CMD build "$p" > "build-$n.txt" 2>&1
  t=$(ls -t ${n}_*.tar.gz 2>/dev/null | head -1)
  if [ -n "$t" ]; then "$R" CMD check --as-cran "$t" > "check-$n.txt" 2>&1
  else echo "NO TARBALL for $n" > "check-$n.txt"; fi
  echo "$n: $(grep -a '^Status' check-$n.txt)"
done

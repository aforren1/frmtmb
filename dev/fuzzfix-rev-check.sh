#!/usr/bin/env bash
# R CMD check --as-cran on core, built with vignettes inside
# dev/fuzzfix-rev-check/ (reviewer copy of dev/fuzzfix-check.sh; the settings of dev/release/run-check.ps1)
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
OUT="$ROOT/dev/fuzzfix-rev-check"
export R_LIBS="C:/Users/adf44/source/r/wt-fuzzfix-rev-lib;C:/Users/adf44/source/r/rellib-r10;C:/Users/adf44/AppData/Local/R/win-library/4.6"
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
export PATH="/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools:/c/Users/adf44/AppData/Roaming/TinyTeX/bin/windows:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
: "${TMP:?TMP is unset}"
mkdir -p "$OUT"; cd "$OUT"
for n in "${@:-frmtmb}"; do
rm -rf "$n.Rcheck" ${n}_*.tar.gz
if [ "$n" = frmtmb ]; then p="$ROOT"; else p="$ROOT/extensions/$n"; fi
"$R" CMD build "$p" > "build-$n.txt" 2>&1
t=$(ls -t ${n}_*.tar.gz 2>/dev/null | head -1)
if [ -n "$t" ]; then "$R" CMD check --as-cran "$t" > "check-$n.txt" 2>&1
else echo "NO TARBALL for $n" > "check-$n.txt"; fi
echo "$n: $(grep -a "^Status" check-$n.txt)"
done

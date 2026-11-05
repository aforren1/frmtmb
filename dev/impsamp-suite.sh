#!/usr/bin/env bash
# One package's whole suite on this lane's library, one R process per
# file, N at a time:  bash dev/impsamp-suite.sh <package> <test dir> [N]
# Logs: dev/impsamp-out/suite-<package>/<file>.txt; RESULT lines on stdout.
pkg="$1"
dir="$2"
N="${3:-6}"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
export R_MAKEVARS_USER="C:/Users/adf44/Documents/.R/Makevars.win"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
out="dev/impsamp-out/suite-$pkg"
mkdir -p "$out"
ls "$dir"/test-*.R | xargs -P "$N" -I{} bash -c '
  f="{}"
  "'"$R"'" dev/impsamp-testfile.R C:/Users/adf44/source/r/wt-impsamp-lib "'"$pkg"'" "$f" \
    > "'"$out"'/$(basename "$f").txt" 2>&1 ||
    echo "FAILED $f exit $?" >> "'"$out"'/$(basename "$f").txt"'
grep -h "^RESULT\|^FAILED" "$out"/*.txt

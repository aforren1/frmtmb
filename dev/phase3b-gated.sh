#!/usr/bin/env bash
# Both extensions' suites with every gated environment variable that
# dev/release/run-gated.ps1 sets, one file per process.
# Usage: bash dev/phase3b-gated.sh <pkg> <base|new>
# Output: dev/phase3b-log/suite-<pkg>-gated-<tag>.txt
cd "$(dirname "$0")/.." || exit 1
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
export FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
export FRMTMB_STAN_CACHE="$(pwd)/dev/stan-cache"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
# "$BASH", not "bash": with Rtools first on PATH, "bash" is Rtools' own
# MSYS bash, and a variable exported here did not reach it.
"$BASH" dev/phase3b-suite.sh "$1" "$2" C:/Users/adf44/source/r/phase3b-lib gated

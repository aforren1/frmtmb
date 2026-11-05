#!/usr/bin/env bash
# Runs one review script with the lane's library stack, in the worktree.
#   bash dev/review-gpgap-run.sh <script> <log> [args...]
set -u
s="$1"; l="$2"; shift 2
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rt=/c/rtools45
export PATH="$PATH:$rt/usr/bin:$rt/x86_64-w64-mingw32.static.posix/bin"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$(pwd)/dev/review-gpgap-cache"
export NOT_CRAN=true
export FRMTMB_BRMS_FIT_TESTS=true
"$RS" "$s" "$@" > "dev/review-gpgap-out/$l" 2>&1
echo "exit $?" >> "dev/review-gpgap-out/$l"

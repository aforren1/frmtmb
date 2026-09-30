#!/usr/bin/env bash
# One test file in one R process: bash dev/formrobust-one.sh <arm> <pkg> <file>
# arm "after" reads the lane library first, "before" reads rellib-r4 only.
# Output: dev/formrobust-log/<arm>--<pkg>--<file>.txt
set -u
arm="$1"; p="$2"; f="$3"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-formrobust
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
export NOT_CRAN=true
if [ "$arm" = "before" ]; then export FORMROBUST_LIB=""; else unset FORMROBUST_LIB; fi
o="$ROOT/dev/formrobust-log/$arm--$p--$(basename "$f" .R).txt"
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$ROOT/dev/formrobust-run-tests.R" "$p" "$f" > "$o" 2>&1
grep -a "^RESULT\|^lib:" "$o"

#!/usr/bin/env bash
# rel084 release review: run one R script from the release tree's root.
#   bash dev/reviews/rel084-review-run.sh <lib r23|r22> <log name> <script> [args]
# Log: dev/reviews/rel084-review-out/<log name> (gitignored).
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export TMP=C:/Users/adf44/AppData/Local/Temp/2
export TEMP=C:/Users/adf44/AppData/Local/Temp/2
export NOT_CRAN=true
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export R_LIBS_USER=C:/Users/adf44/AppData/Local/R/win-library/4.6
export REVLIB=$1; shift
export R_LIBS=C:/Users/adf44/source/r/rellib-$REVLIB
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel084
OUT=$ROOT/dev/reviews/rel084-review-out
log=$1; shift
cd "$ROOT" || exit 1
echo "START $(date '+%F %T') REVLIB=$REVLIB $*" > "$OUT/$log"
"$RS" "$@" >> "$OUT/$log" 2>&1
rc=$?
echo "EXIT $rc $(date '+%F %T')" >> "$OUT/$log"
exit $rc

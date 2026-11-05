#!/usr/bin/env bash
# The rest of fix round 2's runs on the final build, one after another. Both
# library variables take absolute paths: a relative one leaves R on the release
# library in a session, and in the equiv flow fails every pair without a sound.
set -u
cd /c/Users/adf44/source/r/frmtmb-wt-refuse
log=dev/refuse-log/final-steps.txt
: > "$log"
step() { echo "$(date +%H:%M) $1" >> "$log"; }
LIB="C:/Users/adf44/source/r/wt-refuse-lib"
BASE="C:/Users/adf44/source/r/rellib-r16"
step "gated file"
( FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true bash dev/refuse-runtests.sh refuse-gated-out 1 dev/refuse-jobs-l1.txt > dev/refuse-log/gated-l1.txt 2>&1 )
step "build OB"
bash dev/ciharden-openblas.sh 0.3.26 lapack > dev/refuse-log/ob26-build.txt 2>&1
bash dev/ciharden-openblas.sh 0.3.32 lapack > dev/refuse-log/ob32-build.txt 2>&1
step "OB, the gated file"
bash dev/ciharden-scan.sh refuse-ob26l ob0.3.26-lapack-t4 "$LIB" 1 dev/refuse-jobs-l1.txt > dev/refuse-log/ob26.txt 2>&1
bash dev/ciharden-scan.sh refuse-ob32l ob0.3.32-lapack-t4 "$LIB" 1 dev/refuse-jobs-l1.txt > dev/refuse-log/ob32.txt 2>&1
rm -rf dev/ciharden-out
step "equiv"
rm -f dev/refuse-equiv-out/p*.tsv dev/refuse-equiv-out/p*.rds
FRMTMB_EQUIV_LIB="$LIB;$BASE" FORCE=1 bash dev/equiv-run.sh 6 dev/refuse-equiv-out > dev/refuse-log/equiv-final.txt 2>&1
PATH="$PATH:/c/rtools45/usr/bin" "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/refuse-verdict-diff.R > dev/refuse-log/verdict-diff.txt 2>&1
step "sessions"
FRMTMB_RELLIB="$LIB" NOT_CRAN=true bash dev/release/run-sessions.sh refuse-sess 8 > dev/refuse-log/sess-driver.txt 2>&1
step "all done"

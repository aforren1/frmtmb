#!/usr/bin/env bash
# Short check of lane hygiene's fix round: the frmtmb.latent forward
# session once on the lane library, its own R_PARALLEL_PORT, exit code
# kept (the fix round's run died once with no SESSION line).
#   bash dev/reviews/hygiene-short-latent.sh
cd /c/Users/adf44/source/r/frmtmb-wt-hygiene || exit 1
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_RELLIB=C:/Users/adf44/source/r/wt-hygiene-lib NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export R_PARALLEL_PORT=12417
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/reviews/hygiene-short-latent.txt
start=$(date +%s)
"$RS" dev/release/run-session.R frmtmb.latent forward > "$out" 2>&1
rc=$?
echo "EXIT $rc after $(( $(date +%s) - start )) s" >> "$out"

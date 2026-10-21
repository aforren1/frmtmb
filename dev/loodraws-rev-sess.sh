#!/usr/bin/env bash
# Reviewer of lane loodraws: one-session runs of the BASE (main
# checkout's tests, rellib-r10, no FRMTMB_RELLIB), gated and ungated,
# to settle whether 4854e425's sessions leave untagged LEAK lines.
#   bash dev/loodraws-rev-sess.sh <package> <gated|ungated> [order]
# Output: dev/loodraws-rev-out/sess-<package>-<gate>-<order>.txt
set -u
p="$1"; gate="$2"; ord="${3:-forward}"
MAIN=/c/Users/adf44/source/r/frmtmb
WT=/c/Users/adf44/source/r/frmtmb-wt-loodraws
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE="$WT/dev/stan-cache"
export NOT_CRAN=true
unset FRMTMB_RELLIB
if [ "$gate" = gated ]; then export FRMTMB_BRMS_FIT_TESTS=true
else unset FRMTMB_BRMS_FIT_TESTS; fi
if [ "$p" = frmtmb ]; then d="$MAIN/tests/testthat"
else d="$MAIN/extensions/$p/tests/testthat"; fi
"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$MAIN/dev/release/run-session.R" \
  "$p" "$ord" "$d" > "$WT/dev/loodraws-rev-out/sess-$p-$gate-$ord.txt" 2>&1
echo "exit $?"

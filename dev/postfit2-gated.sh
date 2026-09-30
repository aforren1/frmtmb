#!/usr/bin/env bash
# The gated tier of dev/release/run-gated.ps1 (less frmtmb.learn, which
# this lane does not change), on this lane's library, one R process per
# file, 10 at a time. The brms-suite files were run separately
# (dev/postfit2-log/gated-*.txt). Logs: dev/postfit2-gated/.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-postfit2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$WT/dev/postfit2-gated
mkdir -p "$OUT"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
jobs=()
for f in "$WT"/tests/testthat/test-bcm-*.R; do jobs+=("frmtmb $(basename "$f")"); done
for f in test-brms-agreement.R test-brms-likelihood.R test-brms-methods.R \
         test-brms-priors.R test-brms-port.R test-rl-example.R \
         test-drmtmb-agreement.R test-fuzz.R; do
  jobs+=("frmtmb $f")
done
jobs+=("frmtmb.sample test-loo.R" "frmtmb.sample test-sampling-ported.R")
# The gates travel as --gated, which the runner turns into its own
# Sys.setenv(): exported here they did not reach R under xargs, and two
# whole runs skipped every gated block (dev/postfit2-gated-run1/)
printf '%s\n' "${jobs[@]}" | xargs -P 10 -L 1 /usr/bin/bash -c \
  '"'"$R"'" '"$WT"'/dev/postfit2-runtest.R "$0" "$1" --gated > '"$OUT"'/"$0"__"$1".txt 2>&1'
echo "jobs: ${#jobs[@]}"
echo "logs with RESULT: $(grep -l '^RESULT' "$OUT"/*.txt | wc -l)"
grep -h '^RESULT' "$OUT"/*.txt
echo "files with a skip:"
grep -h '^RESULT' "$OUT"/*.txt | grep -v 'skipped=0' || true

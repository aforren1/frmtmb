#!/usr/bin/env bash
# Punch round 1: every extension test file that calls frm_bootstrap(),
# simulate() or refit(), on this lane's core (the extensions themselves
# are the base build's, unchanged by this lane), gated, one R process
# per file, 10 at a time, with frmtmb and the package attached as their
# tests/testthat.R attaches them. Logs: dev/postfit2-p1-ext-<arm>/.
#   bash dev/postfit2-ext-suite.sh <lane|base>
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-postfit2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
ARM=${1:-lane}
OUT=$WT/dev/postfit2-p1-ext-$ARM
EXTRA=""
[ "$ARM" = base ] && EXTRA="--base"
mkdir -p "$OUT"
jobs=()
for p in frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.spline; do
  for f in $(grep -l "frm_bootstrap(\|simulate(\|refit(" \
               "$WT"/extensions/$p/tests/testthat/test-*.R); do
    jobs+=("$p $(basename "$f")")
  done
done
printf '%s\n' "${jobs[@]}" | xargs -P 10 -L 1 /usr/bin/bash -c \
  '"'"$R"'" '"$WT"'/dev/postfit2-runtest.R "$0" "$1" --gated --attach '"$EXTRA"' > '"$OUT"'/"$0"__"$1".txt 2>&1'
echo "jobs: ${#jobs[@]}"
echo "logs with RESULT: $(grep -l '^RESULT' "$OUT"/*.txt | wc -l)"
grep -h '^RESULT' "$OUT"/*.txt | awk '{
  for (i = 1; i <= NF; i++) { split($i, kv, "="); v[kv[1]] += kv[2] }
} END { printf "failed=%d error=%d skipped=%d warning=%d passed=%d\n",
  v["failed"], v["error"], v["skipped"], v["warning"], v["passed"] }'
echo "files with a failure, an error or a skip:"
grep -h '^RESULT' "$OUT"/*.txt | grep -v 'failed=0 error=0 skipped=0' || echo "  none"
echo "files without a RESULT line:"
for f in "$OUT"/*.txt; do grep -q '^RESULT' "$f" || echo "  $(basename "$f")"; done

#!/usr/bin/env bash
# Punch round 1: the whole frmtmb.sample suite, ungated and gated, on
# this lane's library, one R process per file, 10 at a time. Logs:
# dev/postfit2-p1-sample/. The gates travel as --gated (see
# dev/postfit2-runtest.R).
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-postfit2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$WT/dev/postfit2-p1-sample
mkdir -p "$OUT"
jobs=()
for f in "$WT"/extensions/frmtmb.sample/tests/testthat/test-*.R; do
  jobs+=("$(basename "$f")")
done
printf '%s\n' "${jobs[@]}" | xargs -P 10 -L 1 /usr/bin/bash -c \
  '"'"$R"'" '"$WT"'/dev/postfit2-runtest.R frmtmb.sample "$0" --gated > '"$OUT"'/"$0".txt 2>&1'
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

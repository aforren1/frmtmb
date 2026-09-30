#!/usr/bin/env bash
# Whole ungated suites of frmtmb and frmtmb.sample on this lane's
# library, one R process per file, 12 at a time. Each file's log goes to
# dev/postfit2-suite/<pkg>__<file>.txt; the summary counts RESULT lines
# so a file that died without one shows up as missing.
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-postfit2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$WT/dev/postfit2-suite
mkdir -p "$OUT"
export NOT_CRAN=true
unset FRMTMB_BRMS_FIT_TESTS
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
jobs=()
for f in "$WT"/tests/testthat/test-*.R; do jobs+=("frmtmb $(basename "$f")"); done
for f in "$WT"/extensions/frmtmb.sample/tests/testthat/test-*.R; do
  jobs+=("frmtmb.sample $(basename "$f")")
done
# Git's bash by path: with Rtools first on PATH a bare `bash` is Rtools'
# own, which drops exported variables (dev/lane-rules.md, shell traps).
# The runner sets NOT_CRAN itself, so this tier's first run was not
# affected; dev/postfit2-gated.sh's was.
printf '%s\n' "${jobs[@]}" | xargs -P 12 -L 1 /usr/bin/bash -c \
  '"'"$R"'" '"$WT"'/dev/postfit2-runtest.R "$0" "$1" > '"$OUT"'/"$0"__"$1".txt 2>&1'
echo "jobs: ${#jobs[@]}"
echo "logs with RESULT: $(grep -l '^RESULT' "$OUT"/*.txt | wc -l)"
grep -h '^RESULT' "$OUT"/*.txt | awk '{
  for (i = 1; i <= NF; i++) { split($i, kv, "="); v[kv[1]] += kv[2] }
} END { printf "failed=%d error=%d skipped=%d warning=%d passed=%d\n",
  v["failed"], v["error"], v["skipped"], v["warning"], v["passed"] }'
echo "files with failures or errors:"
grep -h '^RESULT' "$OUT"/*.txt | grep -v 'failed=0 error=0' || true
echo "files without a RESULT line:"
for f in "$OUT"/*.txt; do grep -q '^RESULT' "$f" || echo "  $(basename "$f")"; done

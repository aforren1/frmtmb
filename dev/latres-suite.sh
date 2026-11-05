#!/usr/bin/env bash
# One R process per test file, N at a time, against the lane's library.
#   bash dev/latres-suite.sh <tag> <N> [pkg ...]
# "core" names core; an extension by its package name. FRMTMB_BRMS_FIT_TESTS
# comes from the caller. Output: dev/latres-out/<tag>/<pkg>--<file>.txt and
# dev/latres-out/<tag>.summary with one RESULT line per file.
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
tag="$1"; N="$2"; shift 2
out="dev/latres-out/$tag"
mkdir -p "$out"
for p in "$@"; do
  if [ "$p" = core ]; then
    dir=tests/testthat; pk=frmtmb
  else
    dir="extensions/$p/tests/testthat"; pk="$p"
  fi
  for f in "$dir"/test-*.R; do
    while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
    b="$(basename "$f" .R)"
    ( "$RSCRIPT" dev/latres-test1.R "$f" "$pk" > "$out/$pk--$b.txt" 2>&1 ) &
  done
done
wait
grep -h "^RESULT\|^lib:" "$out"/*.txt | grep "^RESULT" > "dev/latres-out/$tag.summary"
n=$(ls "$out"/*.txt | wc -l)
r=$(wc -l < "dev/latres-out/$tag.summary")
echo "ran $r of $n files"
awk '{for(i=3;i<=NF;i++){split($i,kv,"="); s[kv[1]]+=kv[2]}}
     END{printf "tests=%d pass=%d fail=%d err=%d skip=%d warn=%d\n",s["tests"],s["pass"],s["fail"],s["err"],s["skip"],s["warn"]}' \
  "dev/latres-out/$tag.summary"
awk '{for(i=3;i<=NF;i++){split($i,kv,"="); if ((kv[1]=="fail"||kv[1]=="err") && kv[2]>0) {print; next}}}' \
  "dev/latres-out/$tag.summary"
for f in "$out"/*.txt; do
  grep -q "^RESULT" "$f" || echo "NO RESULT LINE: $f"
done

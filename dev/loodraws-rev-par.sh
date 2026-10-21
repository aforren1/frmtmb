#!/usr/bin/env bash
# Reviewer of lane loodraws, fix round 1: core and frmtmb.sample, one
# test file per R process, N at a time, gated, on the reviewer's build
# (wt-loodraws-rev-lib holds the lane's core and frmtmb.sample).
#   bash dev/loodraws-rev-par.sh [N]
# Output: dev/loodraws-rev-out/par/<pkg>--<file>.txt and
# dev/loodraws-rev-out/par-totals.txt
set -u
N="${1:-4}"
WT=/c/Users/adf44/source/r/frmtmb-wt-loodraws
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$WT/dev/loodraws-rev-out/par"
RT=/c/rtools45
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
rm -rf "$OUT"; mkdir -p "$OUT"
jobs_list() {
  for f in "$WT"/tests/testthat/test-*.R; do echo "frmtmb $f"; done
  for f in "$WT"/extensions/frmtmb.sample/tests/testthat/test-*.R; do
    echo "frmtmb.sample $f"
  done
}
jobs_list | while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  ( "$RS" "$WT/dev/loodraws-rev-runtest.R" "$p" "$f" rev \
      > "$OUT/$p--$(basename "$f" .R).txt" 2>&1 ) &
done
wait
{
  for p in frmtmb frmtmb.sample; do
    n=$(ls "$OUT"/$p--*.txt | wc -l)
    r=$(grep -h "^RESULT $p " "$OUT"/$p--*.txt)
    nr=$(echo "$r" | grep -c RESULT)
    echo "$r" | awk -v p=$p -v n=$n -v nr=$nr '{
      for (i = 1; i <= NF; i++) { split($i, kv, "="); s[kv[1]] += kv[2] } }
      END { printf "%s files=%d results=%d pass=%d fail=%d err=%d skip=%d warn=%d\n",
            p, n, nr, s["pass"], s["fail"], s["err"], s["skip"], s["warn"] }'
    grep -L "^RESULT" "$OUT"/$p--*.txt | sed 's/^/NO RESULT: /'
  done
  grep -h "^lib:" "$OUT"/*.txt | sort | uniq -c
  grep -H "^DETAIL \(FAIL\|ERROR\|WARN\)" "$OUT"/*.txt
} > "$WT/dev/loodraws-rev-out/par-totals.txt"
cat "$WT/dev/loodraws-rev-out/par-totals.txt" | head -20

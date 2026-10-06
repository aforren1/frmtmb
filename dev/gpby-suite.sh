#!/usr/bin/env bash
# Every test file of core and the seven extensions, one R process per
# file, N at a time, through dev/gpby-runtest.R.
#   bash dev/gpby-suite.sh <tier: plain|gated> <arm: lane|base> [N]
# Per-file logs in dev/gpby-suite-<tier>-<arm>/, summary in
# dev/gpby-suite-<tier>-<arm>/SUMMARY.txt, generated from the logs.
set -u
tier="$1"; arm="$2"; N="${3:-20}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-gpby
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT="$ROOT/dev/gpby-suite-$tier-$arm"
export NOT_CRAN=true
if [ "$tier" = "gated" ]; then export FRMTMB_BRMS_FIT_TESTS=true; fi
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT"
list="$OUT/jobs.txt"
for f in "$ROOT"/tests/testthat/test-*.R; do echo "frmtmb $f"; done > "$list"
for e in coupling eam latent learn ode sample spline; do
  for f in "$ROOT"/extensions/frmtmb.$e/tests/testthat/test-*.R; do
    echo "frmtmb.$e $f"
  done >> "$list"
done
cd "$ROOT"
while read -r p f; do
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( "$RS" dev/gpby-runtest.R "$p" "$f" $([ "$arm" = base ] && echo base) \
      > "$o" 2>&1 ) &
done < "$list"
wait
{
  echo "tier $tier arm $arm, $(wc -l < "$list") files"
  for p in frmtmb frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn \
           frmtmb.ode frmtmb.sample frmtmb.spline; do
    fs=$(ls "$OUT"/$p--*.txt 2>/dev/null)
    n=0; pa=0; fa=0; er=0; sk=0; wa=0; miss=""
    for o in $fs; do
      n=$((n + 1))
      r=$(grep -a "^RESULT" "$o" | tail -1)
      if [ -z "$r" ]; then miss="$miss $(basename "$o")"; continue; fi
      pa=$((pa + $(echo "$r" | sed 's/.*pass=\([0-9]*\).*/\1/')))
      fa=$((fa + $(echo "$r" | sed 's/.*fail=\([0-9]*\).*/\1/')))
      er=$((er + $(echo "$r" | sed 's/.*error=\([0-9]*\).*/\1/')))
      sk=$((sk + $(echo "$r" | sed 's/.*skip=\([0-9]*\).*/\1/')))
      wa=$((wa + $(echo "$r" | sed 's/.*warn=\([0-9]*\).*/\1/')))
    done
    echo "$p files=$n pass=$pa fail=$fa error=$er skip=$sk warn=$wa"
    [ -n "$miss" ] && echo "  NO RESULT:$miss"
    lib=$(grep -ah "^lib:" $fs 2>/dev/null | sort -u | head -3)
    echo "  $lib"
    for o in $fs; do
      r=$(grep -a "^RESULT" "$o" | tail -1)
      case "$r" in
        *fail=0\ error=0*) ;;
        "") ;;
        *) echo "  BAD $r" ;;
      esac
    done
  done
} > "$OUT/SUMMARY.txt"
cat "$OUT/SUMMARY.txt"

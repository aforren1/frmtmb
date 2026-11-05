#!/usr/bin/env bash
# The files the last edits reach, after the whole suite ran, one per process.
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
out=dev/latres-out/rerun; mkdir -p "$out"
export FRMTMB_BRMS_FIT_TESTS=true
for f in tests/testthat/test-{latres,latres-brms,fcor,importance,message-uniqueness,autocor,autocor-cond,compat,covstruct-seam,guardrail1,registry-helper,tre,extension-home-mixed,compat-extensions,autocor-seam}.R; do
  [ -f "$f" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 1; done
  ( "$RSCRIPT" dev/latres-test1.R "$f" frmtmb > "$out/frmtmb--$(basename "$f" .R).txt" 2>&1 ) &
done
for f in extensions/frmtmb.sample/tests/testthat/test-latres.R extensions/frmtmb.sample/tests/testthat/test-compat-preflight.R extensions/frmtmb.spatial/tests/testthat/test-sar.R; do
  p=$(echo "$f" | cut -d/ -f2)
  while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 1; done
  ( "$RSCRIPT" dev/latres-test1.R "$f" "$p" > "$out/$p--$(basename "$f" .R).txt" 2>&1 ) &
done
wait
grep -h "^RESULT" "$out"/*.txt
grep -h -A1 "^FAIL\|^ERROR" "$out"/*.txt | cut -c1-260

#!/usr/bin/env bash
# Reviewer, lane wt-conditions (recheck 1): spot-check suites on the lane
# build, one file per R process, P at a time: every extension test file
# plus the core files the round's refusals touch.
#   bash dev/conditions-rev-spot.sh <P> <outdir>
P=$1; out=$2
mkdir -p "$out"
date +%H:%M:%S > "$out/.start"
{ ls extensions/*/tests/testthat/test-*.R
  for f in test-conditions test-conditions-census test-interop test-sugar \
           test-nl-lexical test-tabular-inputs test-arg-refusal \
           test-brms-families test-brms-suite-families test-predict-newdata \
           test-brms-methods test-brms-priors test-setprior test-prior-compat \
           test-brms-formula-priors test-priors-bounds-grcov test-ce-bands \
           test-ce-facets test-edgecases test-effects test-v07 test-famgaps \
           test-custom-family test-structure test-aliased-grouping \
           test-dates test-autocor test-ps test-get-prior-route; do
    [ -f tests/testthat/$f.R ] && echo tests/testthat/$f.R
  done
  ls tests/testthat/test-*famlink*.R 2>/dev/null
} | sort -u > "$out/files.txt"
cat "$out/files.txt" | xargs -P "$P" -I{} bash -c '
  f="{}"; log="'"$out"'/$(echo "$f" | tr "/" "_").log"
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" \
    dev/conditions-rev-sweep-onefile.R lane "$f" '"$out"' > "$log" 2>&1
  l=$(grep -h "^BLOCKS" "$log")
  if [ -z "$l" ]; then l="NO RESULT LINE  $f"; fi
  echo "$l"
' > "$out/summary.txt"
date +%H:%M:%S > "$out/.end"

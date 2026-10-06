#!/usr/bin/env bash
# All eight suites again on the release library, one file per process,
# every gate on, with dev/rel068-profile.R as R_PROFILE_USER: it counts
# the fits whose gradient at the optimum is not finite, whether the
# non-finite-gradient warning fired, and whether lane fixes' flat
# warning fired on the same fit (the ordmix final review's rerun list).
#
#   bash dev/rel068-prof-par.sh [N]
# Output: dev/release/prof-files/<pkg>--<file>.txt, the per-process
# lines in dev/rel068-log/prof-cc/<pid>.txt, the summary on stdout.
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
N="${1:-24}"
list=$ROOT/dev/rel068-log/suite.jobs
OUT=$ROOT/dev/release/prof-files
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true
export R_PROFILE_USER="C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel068-profile.R"
export ORDMIX_REV2_LOG="C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel068-log/prof-cc"
: "${TMP:?TMP is unset}"
rm -rf "$OUT" "$ROOT/dev/rel068-log/prof-cc"
mkdir -p "$OUT" "$ROOT/dev/rel068-log/prof-cc"
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/$p--$(basename "$f" .R).txt"
  ( ORDMIX_REV2_FILE="$p/$(basename "$f")" "$RS" \
      "$ROOT/dev/release/run-tests.R" "$p" "$f" > "$o" 2>&1 ) &
done < "$list"
wait
# find, not a glob: 353 long paths pass the Windows command-line limit,
# and cat then read a truncated list
ran=$(find "$OUT" -name "*.txt" -exec grep -h "^RESULT" {} + | grep -c "^RESULT")
echo "PROF ran $ran of $(grep -c '' "$list")"
find "$ROOT/dev/rel068-log/prof-cc" -name "*.txt" -exec cat {} + > "$ROOT/dev/rel068-log/prof-cc.txt"
echo "fits checked: $(grep -c '^CC' "$ROOT/dev/rel068-log/prof-cc.txt")"
echo "non-finite gradient: $(grep -c '^CC.*nonfinite=TRUE' "$ROOT/dev/rel068-log/prof-cc.txt")"
echo "  of them warned: $(grep -c '^CC.*nonfinite=TRUE warned=TRUE' "$ROOT/dev/rel068-log/prof-cc.txt")"
echo "warned without a non-finite gradient: $(grep -c '^CC.*nonfinite=FALSE warned=TRUE' "$ROOT/dev/rel068-log/prof-cc.txt")"
echo "flat warnings: $(grep -c '^FLAT' "$ROOT/dev/rel068-log/prof-cc.txt")"
both=0
for f in "$ROOT"/dev/rel068-log/prof-cc/*.txt; do
  b=$(awk '/^FLAT/{flat=1; next} /^CC/{if (flat && /warned=TRUE/) n++; flat=0} END{print n+0}' "$f")
  both=$((both + b))
done
echo "fits where both warnings fired: $both"
echo "PROF DONE"

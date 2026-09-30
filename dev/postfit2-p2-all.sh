#!/usr/bin/env bash
# Punch round 2: every tier, on this lane's library, one R process per
# file. The gates travel as --gated (see dev/postfit2-runtest.R).
#   bash dev/postfit2-p2-all.sh <tier>
# Tiers (logs under dev/postfit2-p2-<tier>/):
#   target  the ce, postfit and bootstrap files of core and sample, gated
#   core    every core file, ungated
#   gated   core's gated tier (as dev/postfit2-gated.sh)
#   sample  every frmtmb.sample file, gated
#   ext     every file of every other extension, gated, attached
set -u
WT=/c/Users/adf44/source/r/frmtmb-wt-postfit2
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
TIER=$1
OUT=$WT/dev/postfit2-p2-$TIER
mkdir -p "$OUT"
rm -f "$OUT"/*.txt
jobs=()
case $TIER in
  target)
    for f in test-boot.R test-brms-methods.R test-brms-utilities.R \
             test-ce-bands.R test-ce-facets.R test-ce-levels.R \
             test-ce-options.R test-conditional-smooths.R test-nlf.R \
             test-simulate-density.R test-simulate-ergonomics.R \
             test-simulate-newdata.R test-trunc-postfit.R; do
      jobs+=("frmtmb $f --gated")
    done
    for f in test-brms-shapes-draws.R test-conditional-effects-draws.R \
             test-draws-methods.R test-gr-by-draws.R test-postfit-draws.R \
             test-re-formula-draws.R; do
      jobs+=("frmtmb.sample $f --gated")
    done
    P=12 ;;
  core)
    for f in "$WT"/tests/testthat/test-*.R; do
      jobs+=("frmtmb $(basename "$f") --none")
    done
    P=14 ;;
  gated)
    for f in "$WT"/tests/testthat/test-bcm-*.R; do
      jobs+=("frmtmb $(basename "$f") --gated")
    done
    for f in test-brms-agreement.R test-brms-likelihood.R \
             test-brms-methods.R test-brms-priors.R test-brms-port.R \
             test-rl-example.R test-drmtmb-agreement.R test-fuzz.R \
             test-brms-suite-brm.R test-brms-suite-brmsfit-helpers.R \
             test-brms-suite-brmsformula.R test-brms-suite-brmsterms.R \
             test-brms-suite-data-helpers.R test-brms-suite-emmeans.R \
             test-brms-suite-families.R test-brms-suite-methods.R \
             test-brms-suite-priors.R test-brms-suite-standata.R; do
      jobs+=("frmtmb $f --gated")
    done
    P=10 ;;
  sample)
    for f in "$WT"/extensions/frmtmb.sample/tests/testthat/test-*.R; do
      jobs+=("frmtmb.sample $(basename "$f") --gated")
    done
    P=10 ;;
  ext)
    for p in frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn \
             frmtmb.ode frmtmb.spline; do
      for f in "$WT"/extensions/$p/tests/testthat/test-*.R; do
        jobs+=("$p $(basename "$f") --gated")
      done
    done
    P=12 ;;
  *) echo "unknown tier $TIER"; exit 1 ;;
esac
# --attach for every package but core: each package's own
# tests/testthat.R attaches frmtmb and itself
CMD='A=""; [ "$0" != frmtmb ] && A="--attach"; G="$2"; [ "$G" = --none ] && G=""; "'"$R"'" '"$WT"'/dev/postfit2-runtest.R "$0" "$1" $G $A > '"$OUT"'/"$0"__"$1".txt 2>&1'
printf '%s\n' "${jobs[@]}" | xargs -P $P -L 1 /usr/bin/bash -c "$CMD"
echo "tier $TIER, jobs: ${#jobs[@]}"
echo "logs with RESULT: $(grep -l '^RESULT' "$OUT"/*.txt | wc -l)"
grep -h '^RESULT' "$OUT"/*.txt | awk '{
  for (i = 1; i <= NF; i++) { split($i, kv, "="); v[kv[1]] += kv[2] }
} END { printf "failed=%d error=%d skipped=%d warning=%d passed=%d\n",
  v["failed"], v["error"], v["skipped"], v["warning"], v["passed"] }'
echo "files with a failure, an error, a warning or a skip:"
grep -h '^RESULT' "$OUT"/*.txt |
  grep -v 'failed=0 error=0 skipped=0 warning=0' || echo "  none"
echo "files without a RESULT line:"
for f in "$OUT"/*.txt; do
  grep -q '^RESULT' "$f" || echo "  $(basename "$f")"
done

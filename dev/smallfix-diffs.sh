#!/usr/bin/env bash
# One diff per item of lane smallfix, so any subset can be taken:
# dev/smallfix-out/<item>.diff. NEWS.md and frmtmb.latent's NEWS.md are
# in news.diff, because the items' bullets share their sections; the
# findings say which bullet is whose. New files are diffed against
# /dev/null.
#
# Three files are changed by two items since fix round 1: S1's handler
# rewrite (keep_missing()) touches R/priors.R and R/curve.R (S3's) and
# R/predict-brms.R (S2's). Their hunks are split by
# dev/smallfix-split-hunks.R ("split:S1:<file>" and
# "split:other:<file>" below).
#
# Then each item diff is applied alone to a clean archive of HEAD, and
# all of them together, whose result must equal this worktree.
#   bash dev/smallfix-diffs.sh
set -u
cd "$(dirname "$0")/.."
out=dev/smallfix-out
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
one() {
  local name=$1; shift
  : > "$out/$name.diff"
  for f in "$@"; do
    case "$f" in
      split:*)
        local part=${f#split:}; local who=${part%%:*}; local file=${part#*:}
        "$RS" dev/smallfix-split-hunks.R "$file" "$who" >> "$out/$name.diff"
        ;;
      *)
        if git ls-files --error-unmatch "$f" > /dev/null 2>&1; then
          git diff -- "$f" >> "$out/$name.diff"
        else
          git diff --no-index -- /dev/null "$f" >> "$out/$name.diff"
        fi
        ;;
    esac
  done
  echo "$name.diff: $(grep -c '^+++ ' "$out/$name.diff") files," \
    "+$(grep -c '^+[^+]' "$out/$name.diff") -$(grep -c '^-[^-]' "$out/$name.diff")"
}
split3="R/priors.R R/curve.R R/predict-brms.R"
s1r=$(git diff --name-only -- R | grep -v -x -e R/curve-feature.R \
        -e R/priors.R -e R/curve.R -e R/predict-brms.R | tr -d '\r')
one S1-saved-fit-packages \
  $s1r split:S1:R/priors.R split:S1:R/curve.R split:S1:R/predict-brms.R \
  tests/testthat/test-extension-home.R \
  tests/testthat/test-extension-home-mixed.R \
  tests/testthat/helper-absent-package.R \
  dev/rel074-helper-absent-package.R \
  extensions/frmtmb.survival/tests/testthat/helper-absent-package.R \
  extensions/frmtmb.spatial/tests/testthat/helper-absent-package.R \
  .github/workflows/absent-helper-identical.yaml \
  .github/workflows/pkgdown.yaml .github/scripts/saved-fit-check.R \
  dev/smallfix-saved-mixed.sh dev/smallfix-predict-probe.R \
  dev/smallfix-methods-probe.R dev/smallfix-handlers.R \
  dev/smallfix-keep-handlers.R dev/smallfix-head-checks.R
one S2-mixture-mvn-tests \
  extensions/frmtmb.latent/tests/testthat/test-mixture-mvn.R \
  extensions/frmtmb.latent/tests/testthat/test-prior-dpar.R \
  split:other:R/predict-brms.R tests/testthat/test-predict-newdata.R \
  dev/smallfix-mvn-ratios.R dev/smallfix-pd-reached.R \
  dev/smallfix-predict-newdata.R
one S3-refusals \
  split:other:R/priors.R split:other:R/curve.R R/curve-feature.R \
  tests/testthat/test-simplex-pars.R tests/testthat/test-dirichlet-prior.R \
  tests/testthat/test-curve.R \
  dev/smallfix-dirichlet-probe.R dev/smallfix-dirichlet-store.R \
  dev/smallfix-theta-probe.R dev/smallfix-brmsrow.R \
  dev/smallfix-partemplate.R man/frm_curve.Rd man/frm_curve_feature.Rd
one S4-ci-skip-gate \
  .github/scripts/skip-gate.sh .github/workflows/R-CMD-check.yaml \
  .github/workflows/test-coverage.yaml \
  tests/testthat/test-families-extreme.R \
  dev/smallfix-gate-replay.sh dev/smallfix-gate-fixtures.sh \
  dev/smallfix-yaml-parse.R dev/smallfix-ascran-env.R dev/smallfix-gev.R
one news NEWS.md extensions/frmtmb.latent/NEWS.md
one harness dev/smallfix-install.sh dev/smallfix-par.sh dev/smallfix-diffs.sh \
  dev/smallfix-split-hunks.R \
  dev/smallfix-check.sh dev/smallfix-covr.R dev/smallfix-sessions.sh \
  dev/smallfix-findings.md
# every changed or new file is in a diff, and only the split files in two
git status --porcelain --untracked-files=all | awk '{print $2}' |
  grep -v '^dev/smallfix-out/\|^dev/smallfix-ci-log/\|^dev/smallfix-check/' \
  | grep -v '^dev/reviews/' | sort > "$out/changed.lst"
grep -h '^+++ b/' "$out"/*.diff | sed 's#^+++ b/##' | sort > "$out/diffed.lst"
echo "changed but in no diff:"; comm -23 "$out/changed.lst" <(sort -u "$out/diffed.lst")
echo "in two diffs (only these may be: $split3):"; uniq -d "$out/diffed.lst"
# apply: each item alone, then all, to a clean archive of HEAD
tmp=$(mktemp -d)
for it in S1-saved-fit-packages S2-mixture-mvn-tests S3-refusals \
          S4-ci-skip-gate news harness; do
  rm -rf "$tmp/t"; mkdir -p "$tmp/t"
  git archive HEAD | tar -x -C "$tmp/t"
  if (cd "$tmp/t" && git apply --whitespace=nowarn "$OLDPWD/$out/$it.diff" \
        2> "$tmp/err"); then
    echo "applies alone: $it"
  else
    echo "DOES NOT APPLY alone: $it"; head -3 "$tmp/err"
  fi
done
rm -rf "$tmp/t"; mkdir -p "$tmp/t"
git archive HEAD | tar -x -C "$tmp/t"
for it in S1-saved-fit-packages S2-mixture-mvn-tests S3-refusals \
          S4-ci-skip-gate news harness; do
  (cd "$tmp/t" && git apply --whitespace=nowarn "$OLDPWD/$out/$it.diff") ||
    echo "DOES NOT APPLY in sequence: $it"
done
neq=0
while read -r f; do
  if ! cmp -s <(tr -d '\r' < "$f") <(tr -d '\r' < "$tmp/t/$f"); then
    echo "differs after applying all: $f"; neq=$((neq + 1))
  fi
done < "$out/changed.lst"
echo "files that differ after applying all: $neq"
rm -rf "$tmp"

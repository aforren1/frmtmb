#!/usr/bin/env bash
# The test files that use importance =, on a build:
#   bash dev/impsamp-tests-run.sh <label> <lib or "">
# One R process per file, up to 4 at once. Run from the worktree root.
label="$1"
lib="${2:-}"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export NOT_CRAN=true
export FRMTMB_STAN_CACHE="$(cd dev/stan-cache && pwd -W)"
mkdir -p "dev/impsamp-out/tests-$label"
core="test-covstruct-seam.R test-gr-by.R test-grad-verdict.R test-importance.R test-me.R test-ps.R test-rl-example.R test-se-exact-hessian.R test-structure.R"
list=""
for f in $core; do
  list="$list frmtmb|tests/testthat/$f"
done
for f in extensions/frmtmb.learn/tests/testthat/test-factorization.R \
  extensions/frmtmb.learn/tests/testthat/test-rlddm-ndt.R \
  extensions/frmtmb.learn/tests/testthat/test-session.R \
  extensions/frmtmb.learn/tests/testthat/test-surface.R; do
  list="$list frmtmb.learn|$f"
done
for f in extensions/frmtmb.spatial/tests/testthat/test-brms-likelihood-car.R \
  extensions/frmtmb.spatial/tests/testthat/test-car.R; do
  list="$list frmtmb.spatial|$f"
done
list="$list frmtmb.sample|extensions/frmtmb.sample/tests/testthat/test-loo.R"
printf '%s\n' $list | xargs -P 4 -I{} bash -c '
  e="{}"; p="${e%%|*}"; f="${e#*|}"
  "'"$R"'" dev/impsamp-testfile.R "'"$lib"'" "$p" "$f" \
    > "dev/impsamp-out/tests-'"$label"'/$(basename "$f").txt" 2>&1'
grep -h "^RESULT" dev/impsamp-out/tests-"$label"/*.txt

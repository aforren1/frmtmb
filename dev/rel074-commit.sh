#!/usr/bin/env bash
# Commit each lane on its branch, then the tested release tree as one
# merge commit with 19b162b9 and the three lanes as parents, and
# fast-forward main. Modeled on the commit script of the 0.73.0 release
# (the session scratchpad's commit-073.sh). NOT RUN by the consolidation:
# run it only on the user's OK.
#
#   bash dev/rel074-commit.sh
#
# Each lane commits its tree except dev/core-boundary.md, the round's
# shared copy, which only the release tree commits (updated). Main holds
# the same untracked dev/core-boundary.md, which a fast-forward to the
# release commit would overwrite, so git refuses; the script checks
# that main's copy is the round's pre-release text (identical to a
# lane's copy but for the mixture_mvn row the user amended) and moves
# it aside to dev/core-boundary.md.pre074 before the fast-forward.
set -eu
R=/c/Users/adf44/source/r
BASE=19b162b9
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
for l in splitsurv mvnlatent spatial release; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-$l is not on $BASE"; exit 1; }
done
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "main is not on $BASE"; exit 1; }
commit_lane() {
  local n="$1" subject="$2" body="$3"
  git -C "$R/frmtmb-wt-$n" add -A
  git -C "$R/frmtmb-wt-$n" reset -q -- dev/core-boundary.md
  git -C "$R/frmtmb-wt-$n" commit -q -F - <<EOF
$subject

$body

$trailer
EOF
  echo "wt-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
commit_lane splitsurv "frmtmb: curve inference into core, frmtmb.survival, frmtmb.spline retired" \
"frm_curve(), frm_curve_deriv() and frm_curve_feature() move from
frmtmb.spline into core, with the family hook post\$curve_check. The new
extension frmtmb.survival holds cox() and cox_baseline() (from core)
and royston_parmar() and rp_floored() (from frmtmb.spline). New seams
frmtmb_family(simplex_pars =, aterm_data =); a saved fit names the
extension its family needs (R/extension-home.R). frmtmb.spline is
retired. Lane record: dev/splitsurv-findings.md; review:
dev/reviews/2026-10-23-splitsurv.md."
commit_lane mvnlatent "frmtmb.latent: mixture_mvn() moves from core; us_sigma() and us_theta()" \
"mixture_mvn() moves from core into frmtmb.latent, its code unchanged
and its 59 test fits bitwise equal. Core exports us_sigma() and its
inverse us_theta(), the parameterization of an unstructured block, and
documents the mixture component interface. Lane record:
dev/mvnlatent-findings.md; review: dev/reviews/2026-10-23-mvnlatent.md."
commit_lane spatial "frmtmb: the covariance-structure seam, guardrail 1, frmtmb.spatial" \
"frmtmb_register_covstruct() and the covariance structure registry;
R/brms-vocabulary.R with the homes and refusals of every brms name, and
test-guardrail1.R; the new extension frmtmb.spatial holds car() and
spde(); a call around a bar term that names no structure is refused;
bare gr(g) is (1 | g). Lane record: dev/spatial-findings.md; review:
dev/reviews/2026-10-23-spatial.md."

REL="$R/frmtmb-wt-release"
git -C "$REL" add -A
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.74.0: the first moves under the core boundary rule

Merges the three reviewed lanes (splitsurv, mvnlatent, spatial) onto
19b162b9. frmtmb 0.74.0; new frmtmb.survival 0.1.0 and frmtmb.spatial
0.1.0; frmtmb.latent 0.7.0, frmtmb.sample 0.21.1, frmtmb.learn 0.7.3;
frmtmb.spline retired; eam, coupling and ode unchanged. Every package
that moved code needs frmtmb 0.74.0 (measured: dev/rel074-floor.sh).
All nine suites with every gate on: 441 files, 29752 assertions, 0
failures; 54 one-session runs clean; the round's 70 test files clean
under OpenBLAS 0.3.26 and 0.3.32; ported brms suite 607 of 823, as at
0.73.0. R CMD check: frmtmb, latent, survival 1 NOTE (V8); spatial,
sample, learn OK.

Release fix, not a lane change (dev/rel074-release-fix.R):
frmtmb.survival's test-saved-fits.R could not build its "not
installed" library stack under R CMD check (1 ERROR). Both extensions
now copy one helper, helper-absent-package.R, with frmtmb.spatial's
two stacks and probe; frmtmb.survival gains an "installed but
unloadable" arm that runs under R CMD check, and only the "not
installed" arm skips there.

Release fixes from the release review (dev/reviews/2026-10-23-release.md,
FIX; dev/rel074-release-review-fix.R): the README names frmtmb.spatial
and the eight companion packages, and its first example runs; the
README and the brms-migration vignette show how to install every
extension (no meta-package, the user's decision); a cox() fit saved by
frmtmb 0.73.0 is refused by name with a request to refit it; NEWS says
to remove frmtmb.spline; the one-session runner watches the
covariance-structure registry; a CI job keeps the two
helper-absent-package.R copies identical; messages, redirects and
comments. Record: dev/round-20261023.md.

$trailer
EOF
)
rel=$(git -C "$REL" commit-tree "$tree" -p "$BASE" \
  -p wt-splitsurv -p wt-mvnlatent -p wt-spatial -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-release now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
# main's untracked copy of the round document, which the release commit
# now tracks
m="$R/frmtmb/dev/core-boundary.md"
if [ -f "$m" ]; then
  if diff -q <(tr -d '\r' < "$m") \
       <(tr -d '\r' < "$R/frmtmb-wt-splitsurv/dev/core-boundary.md" |
         sed 's/^| `mixture_mvn()` | core | `frmtmb.latent` | an `extension_family_homes` entry |$/| `mixture_mvn()` | core | `frmtmb.latent` | none; no string spelling reaches it, and brms has no such family |/') > /dev/null; then
    mv "$m" "$m.pre074"
    echo "moved main's untracked dev/core-boundary.md to dev/core-boundary.md.pre074"
  else
    echo "main's dev/core-boundary.md is not the round's text; not moving it"; exit 1
  fi
fi
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

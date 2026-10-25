#!/usr/bin/env bash
# Commit each lane on its branch, then the tested release tree as one
# merge commit with 4fd134db and the four lanes as parents, and
# fast-forward main. Modeled on dev/rel074-commit.sh. NOT RUN by the
# consolidation: run it only on the user's OK.
#
#   bash dev/rel075-commit.sh
#
# The lanes are on 18b536fc; the release tree and main are on 4fd134db,
# 18b536fc plus one CI commit (pkgdown.yaml, dev/release/build-docs.R).
# Each lane commits its whole tree (git add -A; logs, outputs and check
# trees are gitignored), the reviewers' dev/review-* scripts with it, as
# at 0.74.0; the release tree did not copy those (dev/rel075-merge.sh),
# so the release commit, whose tree is the release worktree's, does not
# carry them. dev/core-boundary.md is tracked since 0.74.0, so no copy
# is moved aside this time.
set -eu
R=/c/Users/adf44/source/r
BASE=18b536fc
HEADC=4fd134db
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
for l in smallfix bhaz brmsnames sarfcor; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-$l is not on $BASE"; exit 1; }
done
[ "$(git -C "$R/frmtmb-wt-rel075" rev-parse --short=8 HEAD)" = "$HEADC" ] ||
  { echo "wt-rel075 is not on $HEADC"; exit 1; }
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$HEADC" ] ||
  { echo "main is not on $HEADC"; exit 1; }
[ -z "$(git -C "$R/frmtmb" status --porcelain)" ] ||
  { echo "main has uncommitted changes; not fast-forwarding over them"; exit 1; }
commit_lane() {
  local n="$1" subject="$2" body="$3"
  git -C "$R/frmtmb-wt-$n" add -A
  git -C "$R/frmtmb-wt-$n" commit -q -F - <<EOF
$subject

$body

$trailer
EOF
  echo "wt-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
commit_lane smallfix "frmtmb: one refusal names every missing package; four small items" \
"S1: a saved fit that needs extension packages this session cannot load
names every one of them in one refusal (stop_extension_parts());
keep_missing() in the 72 handlers that could swallow the refusal (90
call sites); require_fit_packages() at the head of the methods that
reached one part alone; a third copy of helper-absent-package.R; the
pkgdown job runs .github/scripts/saved-fit-check.R. S2: the
mixture_mvn() test tolerances as ratios, the predict(newdata =)
refusal per model. S3: the Dirichlet refusal at the model (BREAKING),
print(curve, digits =). S4: the CI skip gate for core; the GEV test
reads extraDistr. Lane record: dev/smallfix-findings.md; review:
dev/reviews/2026-10-24-smallfix.md."
commit_lane bhaz "frmtmb.survival: brms's bhaz(); frmtmb_register_aterm(option =)" \
"An addition term may set an argument of the family
(frmtmb_register_aterm(option =)). frmtmb.survival registers bhaz() with
it: bhaz(df =, degree =, intercept =, knots =, Boundary.knots =) sets
cox()'s baseline basis, identical to the cox() spelling, and equal to
brms 2.23.0's basis to 6.7e-16; cox() gains knots and Boundary.knots.
brms-36 (interval-censored cox) in dev/upstream-bugs.md. Lane record:
dev/bhaz-findings.md; review: dev/reviews/2026-10-24-bhaz.md."
commit_lane brmsnames "frmtmb: brms's coefficient names for terms the parser aliased" \
"unsub_specials() undoes the parser's .frm_<name> alias after
splitForm(), so exp(x) is expx everywhere brms has it; the prior table
lists and sorts coefficients as brms does, and set_prior(coef =) takes
both spellings as one slot; fixef() orders nonlinear and multivariate
rows as brms does, and fixef() of draws follows it; the dead alias
code is removed. BREAKING. Lane record: dev/brmsnames-1024-findings.md;
review: dev/reviews/2026-10-24-brmsnames.md."
commit_lane sarfcor "frmtmb: fcor(), the residual correlation seam; sar() in frmtmb.spatial" \
"Part A: frmtmb_register_autocor(), core's fcor(), frmtmb.spatial's sar()
(lag and error), frmtmb_autocor_pars(), the bounds = slot for a class
prior. Part B: frmtmb_covstruct_info(), frmtmb_register_covstruct(
full_rank =) read by the weights check, the autoscale message.
brms-36 to brms-38 of the lane (brms-37 to brms-39 after the
consolidation's renumbering). Lane record: dev/sarfcor-findings.md;
review: dev/reviews/2026-10-24-sarfcor.md."

REL="$R/frmtmb-wt-rel075"
git -C "$REL" add -A
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.75.0: brms names, bhaz(), sar() and fcor(), one refusal

Merges the four reviewed lanes (smallfix, bhaz, brmsnames, sarfcor),
all on 18b536fc, onto 4fd134db. frmtmb 0.75.0 (BREAKING: brms's
coefficient names and orders; bhaz() is frmtmb.survival's; the
Dirichlet refusal at the model); frmtmb.survival 0.2.0,
frmtmb.spatial 0.2.0, frmtmb.sample 0.22.0, frmtmb.latent 0.7.1, each
needing frmtmb 0.75.0 (measured: dev/rel075-floor.sh); eam, coupling,
ode and learn unchanged. New seams: frmtmb_register_aterm(option =),
frmtmb_register_autocor(), frmtmb_autocor_pars(),
frmtmb_covstruct_info(), frmtmb_register_covstruct(full_rank =).

Consolidation (dev/rel075-resolve.R): one refusal names every missing
package across families, covariance structures and residual terms;
frmtmb_covstruct_info() lets the refusal through (keep_missing());
autocor_matrix() and frmtmb_autocor_pars() check the whole fit first;
the upstream entries renumbered (brms-36 bhaz's, brms-37 to brms-39
sarfcor's); new tests of a fit with parts from three absent packages
and of a saved car() and sar() fit.

Release fixes from the release review (dev/reviews/2026-10-24-release.md,
FIX; dev/rel075-release-fix.R): recover_data.frmtmb_fit() has no head
check, which emmeans's try() discarded, so emmeans() and ref_grid() of
a fit with a missing part refuse by the package's name again
(emm_basis.frmtmb_fit()), and the saved-fit test calls both;
frmtmb_covstruct_info("sar") names a residual correlation term; one
prior slot per sar() parameter (no internal row in default_prior(), a
prior on it refused; ar()'s double slot filed); a lagsar prior on an
error SAR names errorsar; frmtmb.survival's README and vignette show
bhaz(); the saved-fit step's comments name both fits;
dev/rel075-floor.sh installs only outside the tree and keeps PATH.

Record: dev/round-20261024.md.

$trailer
EOF
)
rel=$(git -C "$REL" commit-tree "$tree" -p "$HEADC" \
  -p wt-smallfix -p wt-bhaz -p wt-brmsnames -p wt-sarfcor -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel075 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

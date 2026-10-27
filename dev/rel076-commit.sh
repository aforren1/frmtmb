#!/usr/bin/env bash
# Commit each lane on its branch, then the tested release tree as one
# merge commit with 70476dc3 and the three lanes as parents, and
# fast-forward main. Modeled on dev/rel075-commit.sh. NOT RUN by the
# consolidation: run it only on the user's OK.
#
#   bash dev/rel076-commit.sh            # commit
#   bash dev/rel076-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write
#
# The lanes, the release tree and main are all on 70476dc3 (frmtmb
# 0.75.0, pushed). Each lane commits its whole tree (git add -A; logs,
# outputs and check trees are gitignored), the reviewers' dev/review-*
# scripts with it, as at 0.75.0, except the loglik reviewer's three
# output folders of ledger tables (dev/review-loglik-fr1-out-ok/,
# -rec-flip/, -rec-ok/: 340 KB of .tsv and .md that dev/.gitignore does
# not catch, scratch of the review's ledger replay). The release tree
# did not copy the lanes' dev/review-* files (dev/rel076-merge.sh), so
# the release commit, whose tree is the release worktree's, does not
# carry those. It does carry the release review's, as 0.75.0's release
# commit did: dev/reviews/2026-10-27-release.md and the 16
# dev/review-rel076-* scripts, written in the release tree; their logs,
# dev/review-rel076-out/, are gitignored (*-out/). The check below stops
# if they are missing.
set -eu
R=/c/Users/adf44/source/r
BASE=70476dc3
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
excl_loglik=(':(exclude)dev/review-loglik-fr1-out-ok'
             ':(exclude)dev/review-loglik-fr1-out-flip'
             ':(exclude)dev/review-loglik-fr1-rec-flip'
             ':(exclude)dev/review-loglik-fr1-rec-ok')
for l in mocode acprior loglik rel076; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-$l is not on $BASE"; exit 1; }
  [ "$(git -C "$R/frmtmb-wt-$l" branch --show-current)" = "wt-$l" ] ||
    { echo "frmtmb-wt-$l is not on branch wt-$l"; exit 1; }
done
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "main is not on $BASE"; exit 1; }
[ "$(git -C "$R/frmtmb" branch --show-current)" = main ] ||
  { echo "the main checkout is not on branch main"; exit 1; }
[ -z "$(git -C "$R/frmtmb" status --porcelain)" ] ||
  { echo "main has uncommitted changes; not fast-forwarding over them"; exit 1; }
# the release review and its scripts go into the release commit
[ -f "$R/frmtmb-wt-rel076/dev/reviews/2026-10-27-release.md" ] ||
  { echo "no release review in wt-rel076"; exit 1; }
nrev=$(find "$R/frmtmb-wt-rel076/dev" -maxdepth 1 -type f -name 'review-rel076-*' |
  wc -l)
[ "$nrev" -ge 1 ] || { echo "no dev/review-rel076-* scripts"; exit 1; }
echo "release review: dev/reviews/2026-10-27-release.md and $nrev dev/review-rel076-* scripts"
paths_of() {   # <lane>: the pathspec its commit takes
  case "$1" in loglik) echo . "${excl_loglik[@]}";; *) echo .;; esac
}
if $DRY; then
  for l in mocode acprior loglik rel076; do
    # shellcheck disable=SC2046
    n=$(cd "$R/frmtmb-wt-$l" && git add -A -n -- $(paths_of "$l") | wc -l)
    echo "wt-$l: git add -A would take $n paths"
    ( cd "$R/frmtmb-wt-$l" && git add -A -n -- $(paths_of "$l") ) |
      sed "s#^#   #" > "$R/frmtmb-wt-rel076/dev/rel076-log/commit-dry-$l.txt"
  done
  echo "main: fast-forward to the release commit possible:" \
       "$(git -C "$R/frmtmb" merge-base --is-ancestor "$BASE" HEAD && echo yes)"
  echo "DRY RUN DONE (no git write)"
  exit 0
fi
commit_lane() {
  local n="$1" subject="$2" body="$3"
  # shellcheck disable=SC2046
  ( cd "$R/frmtmb-wt-$n" && git add -A -- $(paths_of "$n") )
  git -C "$R/frmtmb-wt-$n" commit -q -F - <<EOF
$subject

$body

$trailer
EOF
  echo "wt-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
commit_lane mocode "frmtmb: mo() codes an integer from its minimum, as brms does" \
"BREAKING since 0.15.0: the simplex of an integer mo() variable is sized
from the minimum, the minimum is stored for new data, and negative
integers are taken. New data that are fractional, missing, a factor or
outside the fitted range are refused, as brms's validate_newdata()
refuses them; conditional_effects() and emmeans hold mo() as brms does;
the standard errors of fitted(), predict(), frm_linpred(),
conditional_effects(), emmeans and marginaleffects include the simplex,
and are NaN at a category no row has; a saved 0.75.0 fit of a mo()
variable that does not start at 0 is refused. brms-40 and brms-41 in
dev/upstream-bugs.md. Lane record: dev/mocode-findings.md; review:
dev/reviews/2026-10-27-mocode.md."
commit_lane acprior "frmtmb: residual correlation priors as brms's" \
"frm_sample()'s defaults are flat on each autocorrelation coefficient,
as brms's are; under a first-order cov = TRUE ar(), ma() or arma() every
prior on class sigma is read at brms's innovation sd; cosy is sampled
on brms's [0, 1] through one logistic chart; a user prior keeps the
bound brms declares for its slot; one prior slot per residual
correlation parameter (the internal theta row is refused), and the prior
tables bound the classes as brms does. New in the sampling API:
autocor_flat_priors(), autocor_bound_chart(), autocor_innov_sd(). Lane
record: dev/acprior-findings.md; review: dev/reviews/2026-10-27-acprior.md."
commit_lane loglik "frmtmb.sample: log_lik() for residual covariance terms" \
"autocor_row_loglik() (sampling API) gives brms's leave-one-out
conditional density of each row given the other rows of its group, read
off the precision matrix; frmtmb_register_autocor(jacobian =) supplies
it, checked against the map at a generic point; frmtmb.spatial's sar()
fills the slot; a dense fallback warns above 1000 rows. frmtmb.sample's
log_lik(), loo(), waic(), psis() and the loo_*() functions take these
models. The 13 bin-2 sar() and fcor() rows of the brms suite are
ported, with the ledger class adapted-input. Lane record:
dev/loglik-findings.md; review: dev/reviews/2026-10-27-loglik.md."

REL="$R/frmtmb-wt-rel076"
git -C "$REL" add -A
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.76.0: mo() coding, residual correlation priors, log_lik()

Merges the three reviewed lanes (mocode, acprior, loglik), all on
70476dc3, onto 70476dc3. frmtmb 0.76.0 (BREAKING: mo() codes an
integer from its minimum; one prior slot per residual correlation
parameter, the prior tables bound those classes as brms does, and a
prior on class sigma under a first-order cov = TRUE ar/ma/arma is read
at the innovation sd); frmtmb.sample 0.23.0 and frmtmb.spatial 0.3.0,
each needing frmtmb 0.76.0 (measured: dev/rel076-floor.sh); survival,
latent, eam, coupling, ode and learn unchanged. New in the sampling
API: autocor_flat_priors(), autocor_bound_chart(), autocor_innov_sd(),
autocor_row_loglik(); new seam slot frmtmb_register_autocor(jacobian =).

Consolidation (dev/rel076-resolve.R): the sampling-API page keeps both
lanes' names; the ported row brmsfit-methods:924 is a divergence by
policy (brms builds the constant mo(Exp) model only without sampling);
new tests of cosy() with mo() and of ar(cov = TRUE) through
frm_sample() and log_lik(). Release fixes (dev/rel076-release-fix.R):
?set_prior says a prior on the sigma intercept under sigma ~ x is still
read at the log marginal sd; frmtmb.sample's NEWS says every prior on
class sigma; the dfbetas() comment; the marginaleffects test's bound is
a ratio to what the simplex adds.

Release review (dev/reviews/2026-10-27-release.md, FIX, no blocker or
major) and its fixes (dev/rel076-release-fix2.R): NEWS lists the refusal
of a two-valued integer mo(); a bound typed at a computed end of a sar()
interval, such as ub = 1, is that end within 64 ulps, and a refusal
prints the ends in 15 digits; ?frm_sample says a prior on the sigma
intercept under sigma ~ x is read at the log marginal sd;
log_lik(newdata =) names posterior_predict(newdata =) only where it
answers; the record's comparison of two tmbstan objectives is corrected
(dev/rtmb-pitfalls.md item 23: tmbstan's process-global pointer).

Record: dev/round-20261027.md.

$trailer
EOF
)
rel=$(git -C "$REL" commit-tree "$tree" -p "$BASE" \
  -p wt-mocode -p wt-acprior -p wt-loglik -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel076 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

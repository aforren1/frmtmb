#!/usr/bin/env bash
# Commit each lane on its branch, then the tested release tree as one
# merge commit with 435e7ce6 and the four lanes as parents, and
# fast-forward main. Modeled on dev/rel076-commit.sh. NOT RUN by the
# consolidation: run it only on the user's OK. It pushes nothing.
#
#   bash dev/rel077-commit.sh            # commit
#   bash dev/rel077-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write
#
# The lanes, the release tree and main are all on 435e7ce6 (frmtmb
# 0.76.0, committed on main, not pushed). Each lane commits its whole
# tree (git add -A; logs, outputs and check trees are gitignored), the
# reviewers' dev/review-* scripts with it, as at 0.76.0; no lane tree
# holds a reviewer output folder that dev/.gitignore lets through (the
# dry run lists every path, and the check below stops on one). The
# release tree did not copy the lanes' dev/review-* files
# (dev/rel077-merge.sh) except the two the backlog cites as the repro of
# the rescor posterior_predict() defect (dev/review-mointerop-rescorpp.R
# and dev/review-mointerop-common.R), so the release commit, whose tree
# is the release worktree's, carries those two and the four reviews in
# dev/reviews/.
set -eu
# git status and similar commands refresh the index's stat data when
# they can take its lock, which is a write; with optional locks off
# they do not, so --dry-run writes nothing (the release review,
# dev/reviews/2026-10-28-release.md, finding 4). git add and git
# commit take the lock regardless.
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=435e7ce6
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
LANES="sigmainn mointerop reslot priordraws"
for l in $LANES rel077; do
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
for l in $LANES; do
  [ -f "$R/frmtmb-wt-rel077/dev/reviews/2026-10-27-$l.md" ] ||
    { echo "no review of $l in wt-rel077"; exit 1; }
done
# a reviewer's output folder that git would take: a directory under
# dev/ named review-*; scripts are files
for l in $LANES; do
  d=$(cd "$R/frmtmb-wt-$l" && git status --porcelain --untracked-files=normal |
        sed -E 's/^.. //' | grep -E '^dev/review-[^/]*/$' || true)
  [ -z "$d" ] || { echo "wt-$l: reviewer output folder caught by git: $d"; exit 1; }
done
# the release review and its fix script go into the release commit; its
# scratch, dev/review-rel077-*-run/, is gitignored (*-run/)
for f in dev/reviews/2026-10-28-release.md dev/rel077-release-fix2.R; do
  [ -f "$R/frmtmb-wt-rel077/$f" ] || { echo "no $f in wt-rel077"; exit 1; }
done
if (cd "$R/frmtmb-wt-rel077" && git add -A -n -- . | grep -q 'dev/review-rel077-'); then
  echo "release review scratch would be committed"; exit 1
fi
echo "release review: dev/reviews/2026-10-28-release.md and dev/rel077-release-fix2.R"
if $DRY; then
  mkdir -p "$R/frmtmb-wt-rel077/dev/rel077-log"
  for l in $LANES rel077; do
    n=$(cd "$R/frmtmb-wt-$l" && git add -A -n -- . | wc -l)
    echo "wt-$l: git add -A would take $n paths"
    ( cd "$R/frmtmb-wt-$l" && git add -A -n -- . ) |
      sed "s#^#   #" > "$R/frmtmb-wt-rel077/dev/rel077-log/commit-dry-$l.txt"
  done
  echo "main: fast-forward to the release commit possible:" \
       "$(git -C "$R/frmtmb" merge-base --is-ancestor "$BASE" HEAD && echo yes)"
  echo "DRY RUN DONE (no git write)"
  exit 0
fi
commit_lane() {
  local n="$1" subject="$2" body="$3"
  ( cd "$R/frmtmb-wt-$n" && git add -A -- . )
  git -C "$R/frmtmb-wt-$n" commit -q -F - <<EOF
$subject

$body

$trailer
EOF
  echo "wt-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
commit_lane sigmainn "frmtmb: sigma is brms's innovation sd under cov = TRUE ar/ma/arma" \
"BREAKING: under the covariance form of ar(), ma() and arma() (first and
higher order), the model's sigma is brms's innovation sd; the likelihood,
the simulator, log_lik() and the residuals read the marginal sd through
autocor_marg_sigma(). The 0.76.0 prior layer that read a sigma prior at
the innovation sd goes (autocor_innov(), innov_logc(), the export
autocor_innov_sd(), prior_logdens(lc =)). conditional_effects(method =
\"posterior_predict\") draws its band at the innovation sd, as brms. A
fit saved by 0.76.0 is refused where an answer reads sigma
(autocor_check_fit(), autocor_check_dpar()); .get_sigma.frmtmb_fit on
insight's internal generic. Ports log_lik:130/133 and
posterior_predict:93/97 (adapted input). Lane record:
dev/sigmainn-findings.md; review: dev/reviews/2026-10-27-sigmainn.md."
commit_lane mointerop "frmtmb: one- and two-value mo(), marginaleffects vcov =" \
"BREAKING: mo() takes a variable with one or two values, as brms does
(Jmo = 0 and 1); the ported block brmsfit-methods:913-949. marginaleffects'
vcov = reaches the fit (a matrix, a function, FALSE, HC0, HC1, stata, a
cluster formula), with the conditional bread cluster_sandwich_full();
gap-aware standard errors through the cumulative simplex weights;
get_predict() passes re_formula, allow_new_levels, dpar and resp and
refuses other packages' arguments by name. Ports the rescor rows
log_lik:107/109 and posterior_predict:69/72. Lane record:
dev/mointerop-findings.md; review: dev/reviews/2026-10-27-mointerop.md."
commit_lane reslot "frmtmb: one prior slot per random-effect and rescor parameter" \
"BREAKING: each random-effect parameter and the residual correlation
have one prior slot, as in brms (cov_class_coefs(), cov_class_refusal());
no bound on a single correlation; the sampling table prints (flat) with
lb/ub for ar, ma, cosy and lagsar; prior = \"flat\" is brms's \"\" in every
slot (kind \"flat\", lkj_dist(flat = TRUE)); the named-list prior
spelling is refused; frm_sample() refuses a route = \"fit\" table
(with_prior_route(), combined_prior_route()); a centered intercept and
Intercept with nlpar have one slot. Ledger priors:5 and :123 pass. Lane
record: dev/reslot-findings.md; review: dev/reviews/2026-10-27-reslot.md."
commit_lane priordraws "frmtmb.sample: sample_prior, prior_draws(), stored warmup" \
"frm_sample(sample_prior = \"no\"/\"yes\"/\"only\") as brms: prior draws
under brms's names, the prior of the whole model through core's new
prior_only frame field in build_objective(); prior_draws(); the warmup
kept once in the stanfit and read with inc_warmup = TRUE; update() of
draws; BREAKING: Evid.Ratio from prior draws, as brms's
evidence_ratio(); the prior overlay of plot() of a hypothesis. Ledger
brmsfit-methods:794, 797 and 798 pass. Lane record:
dev/priordraws-findings.md; review: dev/reviews/2026-10-27-priordraws.md."

REL="$R/frmtmb-wt-rel077"
git -C "$REL" add -A
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.77.0: innovation sigma, one prior slot, mo(), prior draws

Merges the four reviewed lanes (sigmainn, mointerop, reslot,
priordraws), all on 435e7ce6, onto 435e7ce6. frmtmb 0.77.0 (BREAKING:
sigma is brms's innovation sd under cov = TRUE ar/ma/arma; one prior
slot per random-effect and rescor parameter, (flat) as brms's "",
the named-list prior refused; mo() on one or two values; marginaleffects
vcov =); frmtmb.sample 0.24.0, needing frmtmb 0.77.0 (sample_prior,
prior_draws(), stored warmup, Evid.Ratio from prior draws, prior =
"flat" as brms's ""). The other seven extensions keep their versions
(dev/round-20261028.md, "Versions and floors").

Consolidation (dev/rel077-resolve.R): the simplex branch of the prior
resolver keeps both lanes' skips; the bin-2 adapter spec keeps both
lanes' blocks and the tier files are regenerated from it; the dead
innov_sd and overrides reads go; the ledger's 18 rows of brms's own
numeric helpers are divergences by policy. Release fixes
(dev/rel077-release-fix.R): insight::get_sigma() hands a current fit to
insight's default; a coerced fit-route prior table is refused where its
(flat) rows replace brms's defaults; a list of specifications points to
+ and c(); the marginaleffects help on HC1's k and leverages; (flat) on a
mo() simplex is dirichlet(1), as brms's ""; the cross-lane tests.

Release review (dev/reviews/2026-10-28-release.md, FIX, one major)
and its fixes (dev/rel077-release-fix2.R): a row or column selection
and subset() of a prior table keep its route, so frm_sample() still
refuses frm()'s table cut down; a coerced table is marked by its
printed (flat) whatever its source column says; c() takes a prior
table, as brms's does; ?frmtmb-marginaleffects says the sandwich
spellings are refused under a residual correlation term; the commit
script's dry run takes no optional git locks.

Record: dev/round-20261028.md.

$trailer
EOF
)
rel=$(git -C "$REL" commit-tree "$tree" -p "$BASE" \
  -p wt-sigmainn -p wt-mointerop -p wt-reslot -p wt-priordraws -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel077 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

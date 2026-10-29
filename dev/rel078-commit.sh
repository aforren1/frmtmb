#!/usr/bin/env bash
# Commit each lane on its branch, then the tested release tree as one
# merge commit with ac4906bd and the three lanes as parents, and
# fast-forward main. Modeled on dev/rel077-commit.sh. NOT RUN by the
# consolidation: run it only on the user's OK. It pushes nothing.
#
#   bash dev/rel078-commit.sh            # commit
#   bash dev/rel078-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write, no lock taken
#
# The lanes, the release tree and main are all on ac4906bd (frmtmb
# 0.77.0 plus one test fix, pushed). Each lane commits its whole tree
# but the reviewers' output folders that dev/.gitignore lets through (a
# dev/review-*-out*/ folder: equiv's dev/review-equiv-out2/, rerun
# results); logs, .rds files and check trees are gitignored. Every
# `git add` here excludes such a folder, and the script stops if one
# would still be taken. The reviewers' scripts are their evidence and are
# committed: the loose dev/review-meme-* and dev/review-rescorpp-* files
# and equiv's dev/review-equiv-scratch/ (the 0.78.0 release review,
# finding 1). The release tree carries copies of them all
# (dev/rel078-merge.sh), with the four reviews in dev/reviews/.
#
# The dry run takes no lock and writes nothing (the 0.78.0 release
# review, finding 4): `git add -n` takes the index lock even when it adds
# nothing, so the dry run lists the paths read-only instead. Its git
# calls are rev-parse, branch --show-current, status --porcelain,
# ls-files and merge-base --is-ancestor; with GIT_OPTIONAL_LOCKS=0,
# status does not take the lock to refresh the index's stat data, and
# the others only read. git add and git commit, in the commit run, take
# the lock as they must.
set -eu
# a failing git or filter inside a pipeline fails the pipeline, so a
# list built from one is never silently short (the release review's
# short check of fix 2)
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=ac4906bd
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
LANES="equiv meme rescorpp"
# every git add of this script: the whole tree, but no reviewer output
# folder (a directory under dev/ named review-*-out*)
SPEC=(-- . ':(exclude,glob)dev/review-*-out*/**')
OUTDIR='^dev/review-[^/]*-out[^/]*/'
# the paths `git add -A "${SPEC[@]}"` would take in a tree, read-only:
# every modified, deleted or untracked file that is not ignored, from
# git status, checked against git ls-files. It stops the script when a
# listing fails, when the two disagree, or when the list is empty (every
# tree here has work to commit). Called as `lst=$(take_list l) || exit 1`,
# never inside a pipeline, so its refusal reaches the script.
take_list() {
  local a b
  a=$(cd "$R/frmtmb-wt-$1" &&
        git status --porcelain=v1 -z --untracked-files=all "${SPEC[@]}" |
        tr '\0' '\n' | cut -c4- | sort -u) ||
    { echo "wt-$1: git status failed" >&2; return 1; }
  b=$(cd "$R/frmtmb-wt-$1" &&
        git ls-files -z --modified --deleted --others --exclude-standard \
          "${SPEC[@]}" | tr '\0' '\n' | sort -u) ||
    { echo "wt-$1: git ls-files failed" >&2; return 1; }
  [ "$a" = "$b" ] ||
    { echo "wt-$1: git status and git ls-files disagree" >&2; return 1; }
  [ -n "$a" ] ||
    { echo "wt-$1: the list of paths to commit is empty" >&2; return 1; }
  printf '%s\n' "$a"
}
for l in $LANES rel078; do
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
  [ -f "$R/frmtmb-wt-rel078/dev/reviews/2026-10-28-$l.md" ] ||
    { echo "no review of $l in wt-rel078"; exit 1; }
done
for f in dev/round-20261029.md dev/rel078-release-fix.R \
         dev/reviews/2026-10-29-release.md dev/rel078-release-fix2.R \
         extensions/frmtmb.sample/tests/testthat/test-rescor-me.R; do
  [ -f "$R/frmtmb-wt-rel078/$f" ] || { echo "no $f in wt-rel078"; exit 1; }
done
ne=$(ls "$R/frmtmb-wt-rel078/dev/review-equiv-scratch" | wc -l)
[ "$ne" -ge 26 ] ||
  { echo "wt-rel078 holds $ne of equiv's 26 review scripts"; exit 1; }
# no reviewer output folder in any commit
for l in $LANES rel078; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | grep -E "$OUTDIR" || true)
  [ -z "$bad" ] || { echo "wt-$l: reviewer output would be committed:"; echo "$bad"; exit 1; }
done
if $DRY; then
  mkdir -p "$R/frmtmb-wt-rel078/dev/rel078-log"
  for l in $LANES rel078; do
    lst=$(take_list "$l") || exit 1
    printf '%s\n' "$lst" | sed "s#^#   #" \
      > "$R/frmtmb-wt-rel078/dev/rel078-log/commit-dry-$l.txt"
    n=$(grep -c '' "$R/frmtmb-wt-rel078/dev/rel078-log/commit-dry-$l.txt")
    st=$(cd "$R/frmtmb-wt-$l" && git status --porcelain --untracked-files=normal) ||
      { echo "wt-$l: git status failed"; exit 1; }
    left=$(printf '%s\n' "$st" | { grep -E '^\?\? dev/review-[^/]*-out[^/]*/$' || true; } |
             sed -E 's/^\?\? //' | tr '\n' ' ')
    echo "wt-$l: git add would take $n paths; left untracked: ${left:-none}"
  done
  echo "main: fast-forward to the release commit possible:" \
       "$(git -C "$R/frmtmb" merge-base --is-ancestor "$BASE" HEAD && echo yes)"
  echo "DRY RUN DONE (no git write, no lock)"
  exit 0
fi
commit_lane() {
  local n="$1" subject="$2" body="$3"
  ( cd "$R/frmtmb-wt-$n" && git add -A "${SPEC[@]}" )
  git -C "$R/frmtmb-wt-$n" commit -q -F - <<EOF
$subject

$body

$trailer
EOF
  echo "wt-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
commit_lane equiv "dev: a measured equivalence survey against brms 2.23.0" \
"The equivalence survey instrument (dev/equiv-*, dev/equiv-data/): 172
frm_compat() pairs of the 7,387 core pairs, drawn by stratum and
weighted by reach, each fitted by brms 2.23.0 and frmtmb 0.77.0 and
compared at four levels. Headline 83.6 percent, Wald [78.2, 89.0],
reach-weighted agreement with brms; 80.5 with the weak-evidence verdicts
out. Defects D1 to D8 (D1: nl x ordinal drops the body; D4:
frm_curve()/frm_linpred() drop cs() offsets), and 53 compat-row
proposals in dev/equiv-data/compat-proposals.tsv, not applied. No
package code changes. Lane record: dev/equiv-findings.md; review:
dev/reviews/2026-10-28-equiv.md."
commit_lane meme "frmtmb: brms's me() prior classes meanme, sdme, corme" \
"BREAKING: set_prior() takes brms's me() classes meanme, sdme and corme
with brms 2.23.0's meaning on both routes, and default_prior(),
validate_prior() and prior_summary() list brms's rows. A
sampling-route table's group rows print brms's \"\" ((flat)), so an
edited class row reaches the groups (sd, cor, corme, thres(gr =),
delta). frm_sample() samples brms's me() defaults (flat on sdme, lkj(1)
on corme), draws their priors under brms's names, refuses coef with a
bound for the brms classes (three exemptions), reads the
nonlinear-intercept spelling on a ~ 1 as brms's class row, and carries
the grouped corme column. me(x) without sdx answers default_prior().
Lane record: dev/meme-findings.md; review:
dev/reviews/2026-10-28-meme.md."
commit_lane rescorpp "frmtmb.sample: joint rescor draws, brms's multivariate layout" \
"BREAKING: posterior_predict() of a set_rescor(TRUE) model draws the
responses of a row jointly (multivariate normal or t with the shared
nu) through core's new export predict_rescor_draw(); without resp the
draws methods return brms's draws x observations x responses array;
resp takes several responses in the order named; log_lik(resp =) is
the marginal or joint density of the named responses, through core's
rescor_row_loglik(resp =); loo_*() take several responses;
draws_resp_set() refuses a name given twice; predict() refuses to
summarize a category response. Lane record: dev/rescorpp-findings.md;
review: dev/reviews/2026-10-28-rescorpp.md."

REL="$R/frmtmb-wt-rel078"
git -C "$REL" add -A "${SPEC[@]}"
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.78.0: brms's me() prior classes, joint rescor draws

Merges the three reviewed lanes (equiv, meme, rescorpp), all on
ac4906bd, onto ac4906bd. frmtmb 0.78.0 (BREAKING: brms's me() prior
classes meanme, sdme and corme; sampling-route group rows print
brms's (flat); new sampling-API export predict_rescor_draw() and
rescor_row_loglik(resp =)); frmtmb.sample 0.25.0, needing frmtmb
0.78.0 (BREAKING: joint predictive draws of set_rescor(TRUE) models
and brms's multivariate layout; brms's me() defaults and prior draws;
coef with a bound refused). The other seven extensions keep their
versions (dev/round-20261029.md, "Versions and floors"). The
equivalence survey (dev/equiv-*) is measurement only.

Consolidation (dev/rel078-resolve.R): the two NEWS files, the only
files both code lanes change, are rebuilt under one heading each;
brms-43 and brms-44 in dev/upstream-bugs.md; the release harness
reads rellib-r16. Release fixes (dev/rel078-release-fix.R): R1, the
nonlinear-intercept spelling and brms's class row on a ~ 1 resolve
by specificity in either order, as frm() does; R2, the ledger notes of
posterior_predict:69 and :72; the cross-lane test test-rescor-me.R
(mvbf() with rescor and me()).

Release review (dev/reviews/2026-10-29-release.md, FIX, no blocker
or major) and its fixes (dev/rel078-release-fix2.R): prior_summary()
of draws lists both spellings of a nonlinear intercept as written while
the sampler reads the merged slot; bayes_R2() and loo_R2() refuse a
response named twice; the equiv review's evidence scripts are kept;
the commit script's dry run takes no lock; Rtools last on PATH in
dev/release/run-par.sh; a stale comment in the bin-2 adapter spec.

Record: dev/round-20261029.md.

$trailer
EOF
)
rel=$(git -C "$REL" commit-tree "$tree" -p "$BASE" \
  -p wt-equiv -p wt-meme -p wt-rescorpp -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel078 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

#!/usr/bin/env bash
# Commit each lane of the round of 2026-11-06 on its branch, then the
# tested release tree as one merge commit with 3cd2a0c8 and the twelve
# lanes as parents, and fast-forward main. Modeled on dev/rel079-commit.sh.
# NOT RUN by the consolidation: run it only on the user's OK. It pushes
# nothing.
#
#   bash dev/rel080-commit.sh            # commit
#   bash dev/rel080-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write, no lock taken, and
#                                        # each tree's index sha256 stays
#
# The lanes, the release tree and main are all on 3cd2a0c8 (frmtmb 0.79.1
# plus the gp() docs, pushed). Each lane commits its whole tree but the
# machine output that dev/.gitignore lets through: a reviewer's output
# folder (dev/review-*-out*/), a lane's copy of the equivalence survey's
# output (dev/*-equiv-out*/) and a reviewer's trial-merge folder
# (dev/review-*-merge/). Every `git add` here excludes those, and the
# script stops if one would still be taken. The release tree carries
# copies of the lanes' records, reviews and evidence scripts
# (dev/rel080-merge.sh).
#
# main's checkout holds one uncommitted edit,
# .github/workflows/test-coverage.yaml. The ci lane carries that edit and
# more (the skip gate also after a failed run). The script checks that
# the release tree's file already holds main's edit (a three-way merge of
# main's change onto the release file changes nothing), then puts main's
# file back to its HEAD blob and fast-forwards: the edit arrives with the
# release commit. If the check fails, it stops before any write to main.
#
# The dry run takes no lock and writes nothing: its git calls are
# rev-parse, branch --show-current, status --porcelain, ls-files, show,
# merge-file -p (to stdout) and merge-base --is-ancestor, under
# GIT_OPTIONAL_LOCKS=0, and it reads the index's sha256 of every tree
# before and after, and stops if one moved.
set -eu
# a failing git or filter inside a pipeline fails the pipeline, so a
# list built from one is never silently short
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=3cd2a0c8
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
LANES="ci latperf mfxse loglik osa nestquad arstep fitedge bigx nlstart readers superassign"
COV=.github/workflows/test-coverage.yaml
SPEC=(-- . ':(exclude,glob)dev/review-*-out*/**' ':(exclude,glob)dev/*-equiv-out*/**'
      ':(exclude,glob)dev/review-*-merge/**' ':(exclude,glob)dev/*.progress'
      ':(exclude,glob)dev/*.done')
OUTDIR='^dev/(review-[^/]*-out[^/]*|[^/]*-equiv-out[^/]*|review-[^/]*-merge)/'
# the paths `git add -A "${SPEC[@]}"` would take in a tree, read-only:
# every modified, deleted or untracked file that is not ignored, from
# git status, checked against git ls-files. It stops the script when a
# listing fails, when the two disagree, or when the list is empty (every
# tree here has work to commit). Called as `lst=$(take_list l) || exit 1`,
# never inside a pipeline, so its refusal reaches the script.
take_list() {
  local a b
  a=$(cd "$R/frmtmb-wt-$1" &&
        git status --porcelain=v1 -z --untracked-files=all --no-renames "${SPEC[@]}" |
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
# the sha256 of each tree's index, to show the dry run moved none
index_sums() {
  for t in frmtmb $(for l in $LANES rel080; do printf 'frmtmb-wt-%s ' "$l"; done); do
    sha256sum "$(git -C "$R/$t" rev-parse --path-format=absolute --git-path index)" ||
      { echo "$t: no index to read" >&2; return 1; }
  done
}
for l in $LANES; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-r80-$l is not on $BASE"; exit 1; }
  [ "$(git -C "$R/frmtmb-wt-$l" branch --show-current)" = "wt-r80-$l" ] ||
    { echo "frmtmb-wt-$l is not on branch wt-r80-$l"; exit 1; }
done
[ "$(git -C "$R/frmtmb-wt-rel080" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "wt-rel080 is not on $BASE"; exit 1; }
[ "$(git -C "$R/frmtmb-wt-rel080" branch --show-current)" = wt-rel080 ] ||
  { echo "frmtmb-wt-rel080 is not on branch wt-rel080"; exit 1; }
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "main is not on $BASE"; exit 1; }
[ "$(git -C "$R/frmtmb" branch --show-current)" = main ] ||
  { echo "the main checkout is not on branch main"; exit 1; }
# main may hold the coverage edit alone, and only one the release has
mst=$(git -C "$R/frmtmb" status --porcelain) || { echo "main: git status failed"; exit 1; }
case "$mst" in
  "") COVFIX=false ;;
  " M $COV") COVFIX=true ;;
  *) echo "main has uncommitted changes other than $COV; not fast-forwarding"; exit 1 ;;
esac
if $COVFIX; then
  tmp=$(mktemp -d)
  git -C "$R/frmtmb" show "$BASE:$COV" > "$tmp/base"
  tr -d '\r' < "$R/frmtmb/$COV" > "$tmp/main"
  tr -d '\r' < "$R/frmtmb-wt-rel080/$COV" > "$tmp/rel"
  # the release's file is main's edit and more: main's edit applied onto
  # it changes nothing where the two do not overlap, and every line main
  # adds or changes is in it where they do
  git merge-file -p --ours "$tmp/rel" "$tmp/base" "$tmp/main" > "$tmp/merged" ||
    { echo "git merge-file failed on $COV"; rm -rf "$tmp"; exit 1; }
  cmp -s "$tmp/merged" "$tmp/rel" ||
    { echo "main's $COV edit is not in the release tree"; rm -rf "$tmp"; exit 1; }
  diff "$tmp/base" "$tmp/main" | sed -n 's/^> //p' > "$tmp/added" || true
  while IFS= read -r ln; do
    grep -qxF -- "$ln" "$tmp/rel" ||
      { echo "main's line is not in the release's $COV: $ln"; rm -rf "$tmp"; exit 1; }
  done < "$tmp/added"
  # and every line main deletes is gone from the release's file too (one
  # copy fewer than the base holds)
  diff "$tmp/base" "$tmp/main" | sed -n 's/^< //p' > "$tmp/deleted" || true
  kept=$(awk 'FILENAME == ARGV[1] { del[$0]++; next }
              FILENAME == ARGV[2] { nb[$0]++; next }
              { nr[$0]++ }
              END { for (l in del) if (nr[l] > nb[l] - del[l]) print l }' \
           "$tmp/deleted" "$tmp/base" "$tmp/rel")
  [ -z "$kept" ] ||
    { echo "main deletes lines the release's $COV keeps: $kept"; rm -rf "$tmp"; exit 1; }
  rm -rf "$tmp"
  echo "main: its $COV edit is in the release tree; it is put back to HEAD before the fast-forward"
fi
for l in $LANES; do
  [ -f "$R/frmtmb-wt-rel080/dev/reviews/2026-11-06-$l.md" ] ||
    { echo "no review of $l in wt-rel080"; exit 1; }
done
for f in dev/round-20261106.md dev/rel080-merge.sh dev/rel080-resolve.R \
         dev/rel080-release-fix.R dev/rel080-release-fix2.R dev/rel080-release-fix3.R dev/rel080-news.R dev/rel080-backlog.R \
         R/quadrature-nested.R R/joint-sparse.R; do
  [ -f "$R/frmtmb-wt-rel080/$f" ] || { echo "no $f in wt-rel080"; exit 1; }
done
# no machine output folder in any commit
for l in $LANES rel080; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | { grep -E "$OUTDIR" || true; })
  [ -z "$bad" ] || { echo "wt-$l: machine output would be committed:"; echo "$bad"; exit 1; }
done
if $DRY; then
  before=$(index_sums) || exit 1
  mkdir -p "$R/frmtmb-wt-rel080/dev/rel080-log"
  for l in $LANES rel080; do
    lst=$(take_list "$l") || exit 1
    printf '%s\n' "$lst" | sed "s#^#   #" \
      > "$R/frmtmb-wt-rel080/dev/rel080-log/commit-dry-$l.txt"
    n=$(grep -c '' "$R/frmtmb-wt-rel080/dev/rel080-log/commit-dry-$l.txt")
    st=$(cd "$R/frmtmb-wt-$l" && git status --porcelain --untracked-files=normal) ||
      { echo "wt-$l: git status failed"; exit 1; }
    left=$(printf '%s\n' "$st" | { grep -E '^\?\? dev/(review-[^/]*-out[^/]*|[^/]*-equiv-out[^/]*|review-[^/]*-merge)/$' || true; } |
             sed -E 's/^\?\? //' | tr '\n' ' ')
    echo "wt-$l: git add would take $n paths; left untracked: ${left:-none}"
  done
  echo "main: fast-forward to the release commit possible:" \
       "$(git -C "$R/frmtmb" merge-base --is-ancestor "$BASE" HEAD && echo yes)"
  after=$(index_sums) || exit 1
  [ "$before" = "$after" ] ||
    { echo "AN INDEX MOVED IN THE DRY RUN"; exit 1; }
  echo "index sha256 of all 14 trees unchanged"
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
  echo "wt-r80-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
commit_lane ci "CI: the skip gate after a failed run, stage limits, cran-lag" \
"The skip gate runs after a failed check or coverage run; stan-next and
brms-likelihood get step limits; dev/release/cran-lag.R fails closed;
atomic stan-cache writes in four packages' test helpers; two tests that
read platform rounding. Lane record: dev/ci-findings.md; review:
dev/reviews/2026-11-06-ci.md."
commit_lane latperf "frmtmb: a latent residual's standard errors through a sparse factor" \
"BREAKING: fitted(), residuals() and frm_linpred(se.fit = TRUE) on a
latent residual read the joint covariance through the sparse Cholesky
factor of the joint precision and its selected inverse, at any number of
rows (332 s to 0.42 s at 8000 rows). Lane record:
dev/latperf-findings.md; review: dev/reviews/2026-11-06-latperf.md."
commit_lane mfxse "frmtmb: marginaleffects standard errors hold the group-level effects" \
"BREAKING: get_coef() and get_vcov() carry the group-level effects and
the theta a field reads, so the delta method reads the joint covariance;
hypotheses() on the fit tests its coefficients; a user vcov = is carried
into the effects. Lane record: dev/mfxse-findings.md; review:
dev/reviews/2026-11-06-mfxse.md."
commit_lane loglik "frmtmb: logLik() is the data log-likelihood under a prior" \
"BREAKING: logLik() leaves out the log prior; logLik_penalized() gives
the maximized value; df = NA under a density prior, and AIC(), BIC(),
anova() stop there. Lane record: dev/loglik-findings.md; review:
dev/reviews/2026-11-06-loglik.md."
commit_lane osa "frmtmb: OSA residuals read each family's own CDF" \
"BREAKING: residuals(type = \"osa\") defaults to the family's CDF on every
core family that has one, in row blocks of about sqrt(n); weights are
frequencies; TMB-6 filed. Lane record: dev/osa-findings.md; review:
dev/reviews/2026-11-06-osa.md."
commit_lane nestquad "frmtmb: quadrature = TRUE integrates two nested blocks" \
"BREAKING: (1 | school/class) under quadrature = TRUE, by re-anchored
nested Gauss-Hermite rules, with quad_nodes and quad_max_gb; profile
intervals refused on such fits. Lane record: dev/nestquad-findings.md;
review: dev/reviews/2026-11-06-nestquad.md."
commit_lane arstep "frmtmb: an exact gp() beside ar() starts at the model's marginal" \
"The exact gp() start reads the ar() coefficients, closing the step-field
gap; a second descent from the marginal's other maximum off the
gaussian; a constant sigma in the start. Lane record:
dev/arstep-findings.md; review: dev/reviews/2026-11-06-arstep.md."
commit_lane fitedge "frmtmb: gp() in a nonlinear body, ridges, mi() beside a latent" \
"BREAKING: a field in a body's parameter starts at the marginal of the
body made linear; restarts move every body on a ridge; mi() on the
response beside a latent residual; the SE check reads the marginal's
gradient. Lane record: dev/fitedge-findings.md; review:
dev/reviews/2026-11-06-fitedge.md."
commit_lane bigx "frmtmb: a covariate far from zero fits as the covariate centered" \
"BREAKING: the optimizer steps in the coordinates of an ill-conditioned
design from kappa 1e3; the frame's rank test is scaled and centered at
1e-11, and a capped column warns (frmtmb_design_capped). Lane record:
dev/bigx-findings.md; review: dev/reviews/2026-11-06-bigx.md."
commit_lane nlstart "frmtmb: nonlinear fits probe the signs of the body's coefficients" \
"BREAKING: a default nl fit also fits from each sign pattern and returns
the highest point, with a message; nonlinear coefficients step in the
units they move the body by. Lane record: dev/nlstart-findings.md;
review: dev/reviews/2026-11-06-nlstart.md."
commit_lane readers "frmtmb: a new level's mean, se__ on the display scale, nl SEs" \
"BREAKING: fitted() at a new level is the mean over its effects
(sample_new_levels = \"gaussian\"); conditional_effects()'s se__ on the
displayed scale; a nonlinear body has a standard error; lognormal sigma
draws warn frmtmb_mean_not_finite. Lane record: dev/readers-findings.md;
review: dev/reviews/2026-11-06-readers.md."
commit_lane superassign "frmtmb and extensions: no <<- in package code" \
"Every <<- site in core and four extensions writes to an environment it
creates or returns its value; test-no-superassign.R guards it. Results
unchanged. Lane record: dev/superassign-findings.md; review:
dev/reviews/2026-11-06-superassign.md."

REL="$R/frmtmb-wt-rel080"
git -C "$REL" add -A "${SPEC[@]}"
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.80.0: twelve lanes of the round of 2026-11-06

Merges the twelve reviewed lanes (ci, latperf, mfxse, loglik, osa,
nestquad, arstep, fitedge, bigx, nlstart, readers, superassign), all on
3cd2a0c8, onto 3cd2a0c8. frmtmb 0.80.0 (BREAKING: logLik() under a prior,
OSA residuals by CDF, nested quadrature, marginaleffects standard errors
with the group-level effects, a new level's mean, ill-conditioned
designs, nonlinear sign probes). frmtmb.sample 0.27.0 and frmtmb.latent
0.7.3 need frmtmb 0.80.0; frmtmb.ode 0.7.3 and frmtmb.survival 0.2.2.

Consolidation: dev/rel080-merge.sh, dev/rel080-resolve.R (hand
resolutions), dev/rel080-release-fix.R (the brief's merge steps),
dev/rel080-release-fix2.R (conflicts found on the build, test edits),
dev/rel080-release-fix3.R (the release review's fixes),
dev/rel080-release-fix4.R (its short check's text edits),
dev/rel080-news.R, dev/rel080-backlog.R.

Record: dev/round-20261106.md.

$trailer
EOF
)
ps=(-p "$BASE")
for l in $LANES; do ps+=(-p "wt-r80-$l"); done
rel=$(git -C "$REL" commit-tree "$tree" "${ps[@]}" -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel080 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
if $COVFIX; then
  # the release commit carries this edit (checked above)
  git -C "$R/frmtmb" checkout -- "$COV"
fi
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

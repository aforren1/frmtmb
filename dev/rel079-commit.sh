#!/usr/bin/env bash
# Commit each lane on its branch, then the tested release tree as one
# merge commit with d5d0056d and the seven lanes as parents, and
# fast-forward main. Modeled on dev/rel078-commit.sh. NOT RUN by the
# consolidation: run it only on the user's OK. It pushes nothing.
#
#   bash dev/rel079-commit.sh            # commit
#   bash dev/rel079-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write, no lock taken, and
#                                        # each tree's index sha256 stays
#
# The lanes, the release tree and main are all on d5d0056d (frmtmb
# 0.78.0, pushed). Each lane commits its whole tree but the reviewers'
# output folders that dev/.gitignore lets through (a dev/review-*-out*/
# folder); logs, .rds files and check trees are gitignored. Every
# `git add` here excludes such a folder, and the script stops if one would
# still be taken. The release tree carries copies of the lanes' reviews
# and evidence scripts (dev/rel079-merge.sh).
#
# The dry run takes no lock and writes nothing: its git calls are
# rev-parse, branch --show-current, status --porcelain, ls-files and
# merge-base --is-ancestor, under GIT_OPTIONAL_LOCKS=0, and it reads the
# index's sha256 of every tree before and after, and stops if one moved.
set -eu
# a failing git or filter inside a pipeline fails the pipeline, so a
# list built from one is never silently short
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=d5d0056d
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
LANES="refuse impsamp latres gpgap nlfix predread ncgp"
# every git add of this script: the whole tree, but no reviewer output
# folder (a directory under dev/ named review-*-out*)
SPEC=(-- . ':(exclude,glob)dev/review-*-out*/**' ':(exclude,glob)dev/*.progress'
      ':(exclude,glob)dev/*.done')
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
  for t in frmtmb frmtmb-wt-refuse frmtmb-wt-impsamp frmtmb-wt-latres \
           frmtmb-wt-gpgap frmtmb-wt-nlfix frmtmb-wt-predread frmtmb-wt-ncgp \
           frmtmb-wt-rel079; do
    sha256sum "$(git -C "$R/$t" rev-parse --path-format=absolute --git-path index)" ||
      { echo "$t: no index to read" >&2; return 1; }
  done
}
for l in $LANES rel079; do
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
for l in refuse impsamp latres gpgap nlfix predread; do
  [ -f "$R/frmtmb-wt-rel079/dev/reviews/2026-11-02-$l.md" ] ||
    { echo "no review of $l in wt-rel079"; exit 1; }
done
[ -f "$R/frmtmb-wt-rel079/dev/reviews/2026-11-03-ncgp.md" ] ||
  { echo "no review of ncgp in wt-rel079"; exit 1; }
for f in dev/round-20261104.md dev/rel079-release-fix.R dev/rel079-resolve.R \
         dev/rel079-merge.sh R/gp-start.R R/autocor-latent.R; do
  [ -f "$R/frmtmb-wt-rel079/$f" ] || { echo "no $f in wt-rel079"; exit 1; }
done
# no reviewer output folder in any commit
for l in $LANES rel079; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | grep -E "$OUTDIR" || true)
  [ -z "$bad" ] || { echo "wt-$l: reviewer output would be committed:"; echo "$bad"; exit 1; }
done
if $DRY; then
  before=$(index_sums) || exit 1
  mkdir -p "$R/frmtmb-wt-rel079/dev/rel079-log"
  for l in $LANES rel079; do
    lst=$(take_list "$l") || exit 1
    printf '%s\n' "$lst" | sed "s#^#   #" \
      > "$R/frmtmb-wt-rel079/dev/rel079-log/commit-dry-$l.txt"
    n=$(grep -c '' "$R/frmtmb-wt-rel079/dev/rel079-log/commit-dry-$l.txt")
    st=$(cd "$R/frmtmb-wt-$l" && git status --porcelain --untracked-files=normal) ||
      { echo "wt-$l: git status failed"; exit 1; }
    left=$(printf '%s\n' "$st" | { grep -E '^\?\? dev/review-[^/]*-out[^/]*/$' || true; } |
             sed -E 's/^\?\? //' | tr '\n' ' ')
    echo "wt-$l: git add would take $n paths; left untracked: ${left:-none}"
  done
  echo "main: fast-forward to the release commit possible:" \
       "$(git -C "$R/frmtmb" merge-base --is-ancestor "$BASE" HEAD && echo yes)"
  after=$(index_sums) || exit 1
  [ "$before" = "$after" ] ||
    { echo "AN INDEX MOVED IN THE DRY RUN"; exit 1; }
  echo "index sha256 of all 9 trees unchanged"
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
commit_lane refuse "frmtmb: refusals of models brms builds: sds, mm() with |ID|, mi()" \
"BREAKING: brms's class sds on a smooth's standard deviation; a missing
response value no density reads leaves the objective, and the fit is
the complete-case one; mi() on 14 continuous families, integrated out by
Laplace where read, with a warning of class frmtmb_mi_laplace_error;
mm() with |ID| and cov; draws name a value Ymi_<resp>[k]. Lane record:
dev/refuse-findings.md; review: dev/reviews/2026-11-02-refuse.md."
commit_lane impsamp "frmtmb: importance = converges by a proposal that follows theta" \
"The corrected objective's proposal scale moves with the outer parameters,
so a fit settles in fewer rounds: median round count 4 to 1 over 648 fits.
Steps in the fit's own units, a warning where the correction still moves.
BREAKING: a gaussian response returns the Laplace fit, and the stop that
said ROSE is gone. Lane record: dev/impsamp-findings.md; review:
dev/reviews/2026-11-02-impsamp.md."
commit_lane latres "frmtmb: autocor on families with no residual as brms's latent residual" \
"BREAKING: ar(), ma(), arma(), cosy(), unstr() under cov = TRUE on a family
with no residual y - mu carry one latent per row, as brms does; fitted() on
new rows is the mean over the latent; refusals and warnings where the
Laplace approximation fails. brms-45 decided. Lane record:
dev/latres-findings.md; review: dev/reviews/2026-11-02-latres.md."
commit_lane gpgap "frmtmb: the exact gp() L1 gap was the translator's order" \
"brms_gp_order(): the translator reads the positions in brms's first
appearance, which closes the gap of 94 to 1127 nats. The ?frmtmb-gp page,
the nugget warning and a smaller nugget in frm_sample(), whose route ncgp
replaces at the release. Lane record: dev/gpgap-findings.md; review:
dev/reviews/2026-11-02-gpgap.md."
commit_lane nlfix "frmtmb: nonlinear bodies under ordinal, Gamma and fcor" \
"nl = TRUE under the ordinal families, restarts off a saddle, Gamma and
inverse.gaussian starts on a default link, fcor() and autocor on a body,
inv_logit() in a body, beta_binomial's density, the nl x family compat rule.
Lane record: dev/nlfix-findings.md; review: dev/reviews/2026-11-02-nlfix.md."
commit_lane predread "frmtmb: readers of the predictor: cs() layers, me() draws, resp" \
"frm_curve() and frm_linpred() read cs() layers; fitted() on a new me() row is
the mean over the draw; me() draws on new data; empty and several resp in
bayes_R2(), loo_R2(), conditional_effects(), pp_mixture(); a warning of class
frmtmb_interop_se_part. Lane record: dev/predread-findings.md; review:
dev/reviews/2026-11-02-predread.md."
commit_lane ncgp "frmtmb: the exact gp() is non-centered, with a marginal warm start" \
"BREAKING: the exact gp() holds z, brms's zgp, at a nugget of 1e-12 * sdgp^2,
over a farthest-first order, and starts at the model's marginal maximum, or
a working response's off the gaussian. frm_sample() samples the tape frm() fits.
TMB-5 is scoped by gp_tape_scope(). Lane record: dev/ncgp-findings.md; review:
dev/reviews/2026-11-03-ncgp.md."

REL="$R/frmtmb-wt-rel079"
git -C "$REL" add -A "${SPEC[@]}"
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.79.0: the defects the equivalence survey measured

Merges the seven reviewed lanes (refuse, impsamp, latres, gpgap, nlfix,
predread, ncgp), all on d5d0056d, onto d5d0056d. frmtmb 0.79.0 (BREAKING:
non-centered exact gp() with a marginal warm start; latent residuals on
families with no residual; brms's class sds; complete-case mi() rows; nl
under ordinal families; fitted() on a new row is the mean over its draw).
frmtmb.sample 0.26.0 and frmtmb.learn 0.7.4, both needing frmtmb 0.79.0;
frmtmb.eam 0.13.1 and frmtmb.survival 0.2.1 gain compat rows.

Consolidation: dev/rel079-merge.sh takes gpgap against ncgp by whole files,
dev/rel079-resolve.R holds every hand resolution, and
dev/rel079-release-fix.R the reviews' release fixes.

Record: dev/round-20261104.md.

$trailer
EOF
)
ps=(-p "$BASE")
for l in $LANES; do ps+=(-p "wt-$l"); done
rel=$(git -C "$REL" commit-tree "$tree" "${ps[@]}" -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel079 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

#!/usr/bin/env bash
# Commit each lane of the round of 2026-11-09 on its branch, then the
# tested release tree as one merge commit with e0cce906 and the five
# lanes as parents, and fast-forward main. Modeled on dev/rel080-commit.sh.
# NOT RUN by the consolidation: run it only on the user's OK. It pushes
# nothing.
#
#   bash dev/rel081-commit.sh            # commit
#   bash dev/rel081-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write, no lock taken, and
#                                        # each tree's index sha256 stays
#
# The lanes, the release tree and main are all on e0cce906 (frmtmb 0.80.0,
# pushed). Each lane commits its whole tree but the machine output that
# dev/.gitignore lets through: a reviewer's output folder
# (dev/review-*-out*/), a lane's copy of the equivalence survey's output
# (dev/*-equiv-out*/), a reviewer's trial-merge folder
# (dev/review-*-merge/), and lane compat's parallel-race output folders
# (dev/compat-m-race*/, 63 tsv files of a reproduction attempt; its
# findings quote what they showed). Every `git add` here excludes those,
# and the script stops if one would still be taken. The release tree
# carries copies of the lanes' records, reviews and evidence scripts
# (dev/rel081-merge.sh).
#
# main's checkout must be clean: the script stops before any write if it
# is not.
#
# Resuming after a partial run. The lane commits are made one at a time.
# If the script stops after some of them, those lanes' branches are one
# commit past e0cce906, and the HEAD checks stop a rerun. Do not reset
# them. For each lane already committed, check `git -C <tree> log
# --oneline -2` (its commit, then e0cce906) and `git -C <tree> status
# --porcelain` (nothing but the excluded folders). Then run the
# remaining commit_lane calls and the release block by hand, in this
# file's order: the release commit's parents are e0cce906 and the five
# wt-r81-<lane> branches, whichever run made them. If the release
# commit exists and main did not move, run only
# `git -C <main checkout> merge --ff-only <release commit>`.
#
# The dry run takes no lock and writes nothing in git: its git calls are
# rev-parse, branch --show-current, status --porcelain, ls-files and
# merge-base --is-ancestor, under GIT_OPTIONAL_LOCKS=0, and it reads the
# index's sha256 of every tree before and after, and stops if one moved.
# It writes its lists to dev/rel081-log/commit-dry-<tree>.txt.
set -eu
# a failing git or filter inside a pipeline fails the pipeline, so a
# list built from one is never silently short
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=e0cce906
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
LANES="translator compat code8 hollow support"
SPEC=(-- . ':(exclude,glob)dev/review-*-out*/**' ':(exclude,glob)dev/*-equiv-out*/**'
      ':(exclude,glob)dev/review-*-merge/**' ':(exclude,glob)dev/compat-m-race*/**'
      ':(exclude,glob)dev/*.progress' ':(exclude,glob)dev/*.done')
OUTDIR='^dev/(review-[^/]*-out[^/]*|[^/]*-equiv-out[^/]*|review-[^/]*-merge|compat-m-race[^/]*)/'
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
  for t in frmtmb $(for l in $LANES rel081; do printf 'frmtmb-wt-%s ' "$l"; done); do
    sha256sum "$(git -C "$R/$t" rev-parse --path-format=absolute --git-path index)" ||
      { echo "$t: no index to read" >&2; return 1; }
  done
}
# read before the first git call on any tree, so that the dry run's
# "unchanged" covers every read the script makes
before=$(index_sums) || exit 1
for l in $LANES; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-r81-$l is not on $BASE"; exit 1; }
  [ "$(git -C "$R/frmtmb-wt-$l" branch --show-current)" = "wt-r81-$l" ] ||
    { echo "frmtmb-wt-$l is not on branch wt-r81-$l"; exit 1; }
done
[ "$(git -C "$R/frmtmb-wt-rel081" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "wt-rel081 is not on $BASE"; exit 1; }
[ "$(git -C "$R/frmtmb-wt-rel081" branch --show-current)" = wt-rel081 ] ||
  { echo "frmtmb-wt-rel081 is not on branch wt-rel081"; exit 1; }
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "main is not on $BASE"; exit 1; }
[ "$(git -C "$R/frmtmb" branch --show-current)" = main ] ||
  { echo "the main checkout is not on branch main"; exit 1; }
mst=$(git -C "$R/frmtmb" status --porcelain) || { echo "main: git status failed"; exit 1; }
[ -z "$mst" ] || { echo "main has uncommitted changes; not fast-forwarding:"; echo "$mst"; exit 1; }
for l in $LANES; do
  case "$l" in translator) rv=2026-11-09 ;; *) rv=2026-11-10 ;; esac
  [ -f "$R/frmtmb-wt-rel081/dev/reviews/$rv-$l.md" ] ||
    { echo "no review of $l in wt-rel081"; exit 1; }
  [ -f "$R/frmtmb-wt-rel081/dev/$l-findings.md" ] ||
    { echo "no findings of $l in wt-rel081"; exit 1; }
done
for f in dev/round-20261109.md dev/rel081-merge.sh dev/rel081-release-fix.R \
         dev/rel081-release-fix2.R dev/rel081-release-fix3.R \
         dev/rel081-release-fix4.R \
         dev/rel081-news.R dev/rel081-backlog.R \
         dev/rel081-log.md dev/equiv-live-row.R; do
  [ -f "$R/frmtmb-wt-rel081/$f" ] || { echo "no $f in wt-rel081"; exit 1; }
done
# no machine output folder in any commit
for l in $LANES rel081; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | { grep -E "$OUTDIR" || true; })
  [ -z "$bad" ] || { echo "wt-$l: machine output would be committed:"; echo "$bad"; exit 1; }
done
if $DRY; then
  mkdir -p "$R/frmtmb-wt-rel081/dev/rel081-log"
  for l in $LANES rel081; do
    lst=$(take_list "$l") || exit 1
    printf '%s\n' "$lst" | sed "s#^#   #" \
      > "$R/frmtmb-wt-rel081/dev/rel081-log/commit-dry-$l.txt"
    n=$(grep -c '' "$R/frmtmb-wt-rel081/dev/rel081-log/commit-dry-$l.txt")
    st=$(cd "$R/frmtmb-wt-$l" && git status --porcelain --untracked-files=normal) ||
      { echo "wt-$l: git status failed"; exit 1; }
    left=$(printf '%s\n' "$st" | { grep -E '^\?\? dev/(review-[^/]*-out[^/]*|[^/]*-equiv-out[^/]*|review-[^/]*-merge|compat-m-race[^/]*)/$' || true; } |
             sed -E 's/^\?\? //' | tr '\n' ' ')
    echo "wt-$l: git add would take $n paths; left untracked: ${left:-none}"
  done
  echo "main: fast-forward to the release commit possible:" \
       "$(git -C "$R/frmtmb" merge-base --is-ancestor "$BASE" HEAD && echo yes)"
  after=$(index_sums) || exit 1
  [ "$before" = "$after" ] ||
    { echo "AN INDEX MOVED IN THE DRY RUN"; exit 1; }
  echo "index sha256 of all 7 trees unchanged"
  echo "DRY RUN DONE (no git write, no lock)"
  exit 0
fi
# the paths now staged in a tree, as take_list spells them
staged_list() {
  ( cd "$R/frmtmb-wt-$1" &&
      git diff --cached --name-only --no-renames -z | tr '\0' '\n' | sort -u )
}
# git add must stage exactly what take_list listed before it, or the
# commit would hold something the dry run did not show
add_checked() {
  local want got
  want=$(take_list "$1") || exit 1
  ( cd "$R/frmtmb-wt-$1" && git add -A "${SPEC[@]}" )
  got=$(staged_list "$1") || { echo "wt-$1: git diff --cached failed"; exit 1; }
  [ "$want" = "$got" ] ||
    { echo "wt-$1: the staged paths differ from the list; nothing committed here"; exit 1; }
}
commit_lane() {
  local n="$1" subject="$2" body="$3"
  add_checked "$n"
  git -C "$R/frmtmb-wt-$n" commit -q -F - <<EOF
$subject

$body

$trailer
EOF
  echo "wt-r81-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
commit_lane translator "dev: the survey's weak and unmeasured pairs lifted to L1" \
"Translator rules for prefixed smooths, an ordinal mixture's thresholds
and cs() blocks, and brms's tmp_xi; L1 at three fixed points for a
degenerate or collapsed-sd pair, with a per-row audit; the REML marginal
identity; exact gp() pairs scored divergence by decision at L1. dev only.
Lane record: dev/translator-findings.md; review:
dev/reviews/2026-11-09-translator.md."
commit_lane compat "frmtmb: refusals of every addition term no ordinal family takes" \
"frm_compat() derives the ordinal-only refusals from the families'
declarations, for registered terms too; a registered term's key is its
name and argument index. The survey scores against the live
frm_compat() row in a child process, strata frozen. Lane record:
dev/compat-findings.md; review: dev/reviews/2026-11-10-compat.md."
commit_lane code8 "frmtmb: code 8 at the maximum, a boundary stop, a ridge's SE verdict" \
"BREAKING: false convergence (8) at the maximum is the message
frmtmb_converged_at_max; a correlation at its limit is held at its edge
in the convergence check; the SE verdict on a flat ridge no longer
depends on where the fit stops (the 0.80.0 CI red). Lane record:
dev/code8-findings.md; review: dev/reviews/2026-11-10-code8.md."
commit_lane hollow "frmtmb.sample: draws carry \$data; loo() takes the fitted rows" \
"BREAKING: the leave-one-out functions of draws accept newdata that is
exactly the fitted rows; hypothesis() states why Evid.Ratio is NA on a
fit; make_conditions() of draws reads the fit's data. Lane record:
dev/hollow-findings.md; review: dev/reviews/2026-11-10-hollow.md."
commit_lane support "frmtmb: exact densities outside the support for readers" \
"BREAKING: row_lpdf(exact = TRUE), read by log_lik() of draws, gives -Inf
for a GEV row outside its support; the fitted objective keeps its
barrier. A trunc() bound outside the support is no truncation for
exponential(), weibull(), frechet(), shifted_lognormal(), von_mises()
and royston_parmar(). Lane record: dev/support-findings.md; review:
dev/reviews/2026-11-10-support.md."

REL="$R/frmtmb-wt-rel081"
add_checked rel081
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.81.0: five lanes of the round of 2026-11-09

Merges the five reviewed lanes (translator, compat, code8, hollow,
support), all on e0cce906, onto e0cce906. frmtmb 0.81.0 (BREAKING:
code 8 at the maximum is a message, a flat ridge's SE verdict, exact
log_lik() of draws outside a support, trunc() bounds outside a support).
frmtmb.sample 0.28.0 needs frmtmb 0.81.0; frmtmb.survival 0.2.3.

Consolidation: dev/rel081-merge.sh, dev/rel081-release-fix.R (the
brief's merge steps), dev/rel081-release-fix2.R (the upstream bug
number), dev/rel081-release-fix3.R (the release review's five
minors), dev/rel081-release-fix4.R (its short check),
dev/rel081-news.R, dev/rel081-backlog.R.

Record: dev/round-20261109.md.

$trailer
EOF
)
ps=(-p "$BASE")
for l in $LANES; do ps+=(-p "wt-r81-$l"); done
rel=$(git -C "$REL" commit-tree "$tree" "${ps[@]}" -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel081 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

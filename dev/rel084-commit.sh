#!/usr/bin/env bash
# Commit each lane of the round of 2026-11-16 on its branch (the eleven
# code lanes and the three measurement lanes, whose commits hold dev/
# files only), then the tested release tree as one merge commit with
# b255350e and the fourteen lanes as parents, and fast-forward main.
# Modeled on dev/rel083-commit.sh.
# NOT RUN by the consolidation: run it only on the user's OK. It pushes
# nothing.
#
#   bash dev/rel084-commit.sh            # commit
#   bash dev/rel084-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write, no lock taken, and
#                                        # each tree's index sha256 stays
#
# The lanes, the release tree and main are all on b255350e (frmtmb 0.83.0
# plus the pkgcheck CI fix). Each lane commits its whole tree but the
# machine output that dev/.gitignore lets through: a reviewer's output
# folder (dev/review-*-out*/), a lane's copy of the equivalence survey's
# output (dev/*-equiv-out*/), a reviewer's trial-merge folder
# (dev/review-*-merge/), lane compat's race folders of 0.81.0
# (dev/compat-m-race*/), *.progress and *.done, the Python byte-code
# cache a reviewer's script can leave in .github/scripts/__pycache__/, and
# this round's reviewers' copies of base files
# (dev/reviews/sparsese-review-base/,
# dev/reviews/ordstarts-build-review-base/). Every `git add` here excludes
# those, and the script stops if one would still be taken. A measurement
# lane (samplers, codesize, cran) may take dev/ paths only. The release
# tree carries copies of the lanes' records, reviews and evidence scripts
# (dev/rel084-merge.sh).
#
# main's checkout must be clean: the script stops before any write if it
# is not.
#
# Resuming after a partial run. The lane commits are made one at a time.
# If the script stops after some of them, those lanes' branches are one
# commit past b255350e, and the HEAD checks stop a rerun. Do not reset
# them. For each lane already committed, check `git -C <tree> log
# --oneline -2` (its commit, then b255350e) and that bash
# dev/rel084-list.sh <tree> lists nothing but the excluded paths. Then
# run the remaining commit_lane calls and the release block by hand, in
# this file's order: the release commit's parents are b255350e and the
# fourteen wt-r84-<lane> branches, whichever run made them. If the
# release commit exists and main did not move, run only
# `git -C <main checkout> merge --ff-only <release commit>`.
#
# The dry run takes no lock and writes nothing in git: its git calls are
# rev-parse, branch --show-current, ls-files, hash-object (no -w), show,
# merge-file on scratch copies and merge-base --is-ancestor, under
# GIT_OPTIONAL_LOCKS=0; never git status or git diff, which refresh an
# index's stat data. It reads the index's sha256 of every tree before and
# after, and stops if one moved. It writes its lists to
# dev/rel084-log/commit-dry-<tree>.txt.
set -eu
# a failing git or filter inside a pipeline fails the pipeline, so a
# list built from one is never silently short
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=b255350e
RELBASE=b255350e
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
CODE_LANES="sparsese gpselimit slopetrig alapse nlband quadesc nbrobust squant drawdir ordstarts cdfname"
MEASURE_LANES="samplers codesize cran"
LANES="$CODE_LANES $MEASURE_LANES"
SPEC=(-- . ':(exclude,glob)dev/review-*-out*/**' ':(exclude,glob)dev/*-equiv-out*/**'
      ':(exclude,glob)dev/review-*-merge/**' ':(exclude,glob)dev/compat-m-race*/**'
      ':(exclude,glob)dev/*.progress' ':(exclude,glob)dev/*.done'
      ':(exclude,glob).github/scripts/__pycache__/**'
      ':(exclude,glob)dev/reviews/*-review-base/**')
OUTDIR='^(dev/(review-[^/]*-out[^/]*/|[^/]*-equiv-out[^/]*/|review-[^/]*-merge/|compat-m-race[^/]*/|reviews/[^/]*-review-base/)|\.github/scripts/__pycache__/)'
# The paths `git add -A "${SPEC[@]}"` would take in a tree, read-only,
# without git status or git diff (both refresh the index's stat data,
# whatever GIT_OPTIONAL_LOCKS says). First list: git ls-files names the
# untracked files that are not ignored, the deleted files, and the
# tracked files whose stat data moved; of those, a tracked file is kept
# only when its content (git hash-object --stdin-paths, no -w; it applies
# the path's CRLF filter) differs from the blob git ls-files -s holds.
# Second list, made independently: the content of EVERY tracked file
# against its blob, whatever its stat data says. The two must agree on the
# tracked files, or the script stops; it stops too when a listing fails or
# when the list is empty (every tree here has work to commit). Called as
# `lst=$(take_list l) || exit 1`, never inside a pipeline, so its refusal
# reaches the script.
# [dir] (default .): git -C, so that main is read without a shell cd into it
blob_moved() {  # paths on stdin; prints those whose content left the index
  local p g=(git -C "${1:-.}")
  p=$(cat)
  [ -n "$p" ] || return 0
  paste -d '\t' <(printf '%s\n' "$p") <(printf '%s\n' "$p" | "${g[@]}" hash-object --stdin-paths) |
    awk -F'\t' 'NR == FNR { split($1, m, " "); blob[$2] = m[2]; next }
                !($1 in blob) || blob[$1] != $2 { print $1 }' \
      <("${g[@]}" ls-files -s -z | tr '\0' '\n') -
}
take_list() {
  local d a del mod moved all_moved
  d="$R/frmtmb-wt-$1"
  a=$(cd "$d" && git ls-files -z --others --exclude-standard "${SPEC[@]}" | tr '\0' '\n' | sed '/^$/d') ||
    { echo "wt-$1: git ls-files --others failed" >&2; return 1; }
  del=$(cd "$d" && git ls-files -z --deleted "${SPEC[@]}" | tr '\0' '\n' | sed '/^$/d' | sort -u) ||
    { echo "wt-$1: git ls-files --deleted failed" >&2; return 1; }
  mod=$(cd "$d" && git ls-files -z --modified "${SPEC[@]}" | tr '\0' '\n' | sed '/^$/d' | sort -u) ||
    { echo "wt-$1: git ls-files --modified failed" >&2; return 1; }
  [ -z "$del" ] || mod=$(comm -23 <(printf '%s\n' "$mod") <(printf '%s\n' "$del"))
  moved=$(cd "$d" && printf '%s' "$mod" | blob_moved | sort -u) ||
    { echo "wt-$1: git hash-object or git ls-files -s failed" >&2; return 1; }
  all_moved=$(cd "$d" && git ls-files -z "${SPEC[@]}" | tr '\0' '\n' | sed '/^$/d' |
                { if [ -n "$del" ]; then grep -vxF -f <(printf '%s\n' "$del"); else cat; fi; } |
                blob_moved | sort -u) ||
    { echo "wt-$1: the full content compare failed" >&2; return 1; }
  [ "$moved" = "$all_moved" ] ||
    { echo "wt-$1: the stat listing and the content compare disagree" >&2; return 1; }
  a=$(printf '%s\n%s\n%s\n' "$a" "$del" "$moved" | sed '/^$/d' | sort -u)
  [ -n "$a" ] ||
    { echo "wt-$1: the list of paths to commit is empty" >&2; return 1; }
  printf '%s\n' "$a"
}
# the sha256 of each tree's index, to show the dry run moved none
index_sums() {
  for t in frmtmb $(for l in $LANES rel084; do printf 'frmtmb-wt-%s ' "$l"; done); do
    sha256sum "$(git -C "$R/$t" rev-parse --path-format=absolute --git-path index)" ||
      { echo "$t: no index to read" >&2; return 1; }
  done
}
# read before the first git call on any tree, so that the dry run's
# "unchanged" covers every read the script makes
before=$(index_sums) || exit 1
for l in $LANES; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-r84-$l is not on $BASE"; exit 1; }
  [ "$(git -C "$R/frmtmb-wt-$l" branch --show-current)" = "wt-r84-$l" ] ||
    { echo "frmtmb-wt-$l is not on branch wt-r84-$l"; exit 1; }
done
[ "$(git -C "$R/frmtmb-wt-rel084" rev-parse --short=8 HEAD)" = "$RELBASE" ] ||
  { echo "wt-rel084 is not on $RELBASE"; exit 1; }
[ "$(git -C "$R/frmtmb-wt-rel084" branch --show-current)" = wt-rel084 ] ||
  { echo "frmtmb-wt-rel084 is not on branch wt-rel084"; exit 1; }
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$RELBASE" ] ||
  { echo "main is not on $RELBASE"; exit 1; }
[ "$(git -C "$R/frmtmb" branch --show-current)" = main ] ||
  { echo "the main checkout is not on branch main"; exit 1; }
# main must be clean, read as take_list reads a lane (no git status)
M="$R/frmtmb"
mst=$({ git -C "$M" ls-files -z --others --exclude-standard | tr '\0' '\n'
        git -C "$M" ls-files -z --deleted | tr '\0' '\n'
        git -C "$M" ls-files -z --modified | tr '\0' '\n' | sed '/^$/d' | blob_moved "$M"
      } | sed '/^$/d' | sort -u) || { echo "main: git ls-files failed"; exit 1; }
[ -z "$mst" ] || { echo "main has uncommitted changes; not fast-forwarding:"; echo "$mst"; exit 1; }
for l in $CODE_LANES; do
  [ -f "$R/frmtmb-wt-rel084/dev/reviews/2026-11-17-$l.md" ] ||
    { echo "no review of $l in wt-rel084"; exit 1; }
done
for l in $LANES; do
  [ -f "$R/frmtmb-wt-rel084/dev/$l-findings.md" ] ||
    { echo "no findings of $l in wt-rel084"; exit 1; }
done
for f in dev/round-20261116.md dev/rel084-merge.sh dev/rel084-list.sh \
         dev/rel084-release-fix.sh dev/rel084-news.R dev/rel084-log.md \
         dev/cdfname-rename.R; do
  [ -f "$R/frmtmb-wt-rel084/$f" ] || { echo "no $f in wt-rel084"; exit 1; }
done
# no machine output in any commit, and dev/ paths only from a
# measurement lane
for l in $LANES rel084; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | { grep -E "$OUTDIR" || true; })
  [ -z "$bad" ] || { echo "wt-$l: machine output would be committed:"; echo "$bad"; exit 1; }
done
for l in $MEASURE_LANES; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | { grep -v '^dev/' || true; })
  [ -z "$bad" ] || { echo "wt-$l: a measurement lane would commit outside dev/:"; echo "$bad"; exit 1; }
done
# The drift guard (0.82.0's release review, major 2). A lane that moved
# after the merge would be committed on its branch and made a parent of a
# release commit whose tree lacks the move; git would then count the
# change as merged, and no later merge would bring it in. So each lane
# must take the number of paths the merge took (dev/rel084-merge.sh, its
# taken list), and every path that one lane alone takes must equal the
# release tree's copy, line endings aside, but the paths the consolidation
# edited after the merge (CONS_EDITS, listed in dev/round-20261116.md).
declare -A MERGED=([sparsese]=48 [gpselimit]=56 [slopetrig]=41 [alapse]=60
  [nlband]=65 [quadesc]=37 [nbrobust]=40 [squant]=65 [drawdir]=47
  [ordstarts]=103 [cdfname]=71 [samplers]=29 [codesize]=15 [cran]=11)
# extensions/*/NEWS.md of eam, survival, coupling: rebuilt
# (dev/rel084-news.R); R/ald-fixed.R, R/ald-laplace.R: merge step 5;
# man/*.Rd: roxygen output of the merged tree; R/families-quantile.R and
# man/quantile_elf.Rd: the release review's F2 (cite the tracked scripts,
# not their ignored logs)
CONS_EDITS=" extensions/frmtmb.eam/NEWS.md extensions/frmtmb.survival/NEWS.md
  extensions/frmtmb.coupling/NEWS.md R/ald-fixed.R R/ald-laplace.R
  man/frmtmb-families.Rd R/families-quantile.R man/quantile_elf.Rd "
# lane:path pairs whose lane lines the consolidation rewrote: the order
# vector of se_lost_clauses() (merge step 3), the pole-path filter of
# nl_fit_probed() (merge step 2), the gp_pred() rows of frm_compat()
# (merge step 6), and the asym_laplace help section (merge step 5)
MULTI_EXEMPT=" alapse:R/se-check.R slopetrig:R/se-check.R
  nlband:R/par-template.R gpselimit:R/par-template.R
  squant:R/compat.R alapse:R/compat.R alapse:R/families.R "
merges_clean() {  # lane path tmpdir: the lane's copy merges into the release copy unchanged
  local rc=0
  git -C "$R/frmtmb-wt-rel084" show "$BASE:$2" 2>/dev/null | tr -d '\r' > "$3/base" || : > "$3/base"
  tr -d '\r' < "$R/frmtmb-wt-rel084/$2" > "$3/rel"
  tr -d '\r' < "$R/frmtmb-wt-$1/$2" > "$3/lane"
  git merge-file -p "$3/rel" "$3/base" "$3/lane" > "$3/out" 2>/dev/null || rc=$?
  [ "$rc" = 0 ] && cmp -s "$3/out" "$3/rel"
}
# When a lane's copy does not merge into the release copy by itself (a
# hunk of another lane sits next to its own, which git merge-file calls
# a conflict), replay the merge as dev/rel084-merge.sh made it: from
# b255350e's blob, every lane that takes the path, in merge order. The
# release copy must be what the replay gives.
replays_clean() {  # path lanes tmpdir
  local q rc
  git -C "$R/frmtmb-wt-rel084" show "$BASE:$1" 2>/dev/null | tr -d '\r' > "$3/base" || : > "$3/base"
  cp "$3/base" "$3/cur"
  for q in $LANES; do
    case " $2 " in *" $q "*) ;; *) continue ;; esac
    tr -d '\r' < "$R/frmtmb-wt-$q/$1" > "$3/lane"
    rc=0; git merge-file -q "$3/cur" "$3/base" "$3/lane" 2>/dev/null || rc=$?
    [ "$rc" = 0 ] || return 1
  done
  cmp -s "$3/cur" <(tr -d '\r' < "$R/frmtmb-wt-rel084/$1")
}
drift_check() {
  local l lst n all="" single p bad=0 nsingle=0 tmp
  for l in $LANES; do
    lst=$(take_list "$l") || return 1
    n=$(printf '%s\n' "$lst" | grep -c .)
    [ "$n" = "${MERGED[$l]:-none}" ] ||
      { echo "wt-$l: takes $n paths, the merge took ${MERGED[$l]:-none}" >&2; return 1; }
    all="$all$(printf '%s\n' "$lst" | sed "s#^#$l #")"$'\n'
  done
  tmp=$(mktemp -d) || return 1
  single=$(printf '%s' "$all" |
             awk 'NF { c[$2]++; w[$2] = $1 } END { for (p in c) if (c[p] == 1) print w[p], p }' |
             sort)
  while read -r l p; do
    [ -n "$p" ] || continue
    nsingle=$((nsingle + 1))
    case "$CONS_EDITS" in *[[:space:]]"$p"[[:space:]]*) continue ;; esac
    if [ ! -e "$R/frmtmb-wt-rel084/$p" ] ||
         ! cmp -s <(tr -d '\r' < "$R/frmtmb-wt-$l/$p") <(tr -d '\r' < "$R/frmtmb-wt-rel084/$p"); then
      echo "DRIFT wt-$l $p: differs from the release tree" >&2
      bad=$((bad + 1))
    fi
  done <<< "$single"
  [ "$bad" = 0 ] || { rm -rf "$tmp"; echo "$bad single-lane paths drifted; re-merge them" >&2; return 1; }
  # A path several lanes take: each lane's copy, merged three-way into the
  # release copy over b255350e's blob, must give the release copy back, or
  # the lane holds a change the release tree lacks. NEWS files, the
  # backlog and dev/upstream-bugs.md were merged with --union and rebuilt,
  # and the consolidation rewrote lanes' lines in the pairs of MULTI_EXEMPT.
  local multi nmulti=0 nexempt=0 nreplay=0 takers
  multi=$(printf '%s' "$all" |
            awk 'NF { c[$2]++; w[$2] = w[$2] " " $1 } END { for (p in c) if (c[p] > 1) { n = split(w[p], l, " "); for (i = 1; i <= n; i++) print l[i], p } }' |
            sort)
  while read -r l p; do
    [ -n "$p" ] || continue
    case "$p" in NEWS.md|extensions/*/NEWS.md|dev/test-backlog.md|dev/upstream-bugs.md|man/*.Rd|NAMESPACE) nexempt=$((nexempt + 1)); continue ;; esac
    case "$MULTI_EXEMPT" in *[[:space:]]"$l:$p"[[:space:]]*) nexempt=$((nexempt + 1)); continue ;; esac
    nmulti=$((nmulti + 1))
    merges_clean "$l" "$p" "$tmp" || {
      takers=$(printf '%s\n' "$multi" | awk -v p="$p" '$2 == p { print $1 }' | tr '\n' ' ')
      if replays_clean "$p" "$takers" "$tmp"; then nreplay=$((nreplay + 1))
      else echo "DRIFT wt-$l $p: its copy does not merge into the release copy unchanged" >&2; bad=$((bad + 1)); fi; }
  done <<< "$multi"
  rm -rf "$tmp"
  [ "$bad" = 0 ] || { echo "$bad multi-lane pairs drifted; re-merge them" >&2; return 1; }
  echo "drift guard: 14 lane counts as merged; $nsingle single-lane paths equal the release tree (consolidation edits aside); $nmulti multi-lane lane-path pairs merge into the release copy unchanged, $nreplay of them by a replay of the merge ($nexempt exempt)"
}
drift_check || exit 1
if $DRY; then
  mkdir -p "$R/frmtmb-wt-rel084/dev/rel084-log"
  for l in $LANES rel084; do
    lst=$(take_list "$l") || exit 1
    printf '%s\n' "$lst" | sed "s#^#   #" \
      > "$R/frmtmb-wt-rel084/dev/rel084-log/commit-dry-$l.txt"
    n=$(grep -c '' "$R/frmtmb-wt-rel084/dev/rel084-log/commit-dry-$l.txt")
    left=$(cd "$R/frmtmb-wt-$l" && git ls-files -z --others --exclude-standard | tr '\0' '\n' |
             { grep -E "$OUTDIR" || true; } | tr '\n' ' ')
    echo "wt-$l: git add would take $n paths; left untracked: ${left:-none}"
  done
  echo "main: fast-forward to the release commit possible:" \
       "$(git -C "$R/frmtmb" merge-base --is-ancestor "$RELBASE" HEAD && echo yes)"
  after=$(index_sums) || exit 1
  [ "$before" = "$after" ] ||
    { echo "AN INDEX MOVED IN THE DRY RUN"; exit 1; }
  echo "index sha256 of all 16 trees unchanged"
  echo "DRY RUN DONE (no git write, no lock)"
  exit 0
fi
# the paths now staged in a tree, as take_list spells them (plumbing:
# reads the index against HEAD, refreshes nothing)
staged_list() {
  ( cd "$R/frmtmb-wt-$1" &&
      git diff-index --cached --name-only --no-renames -z HEAD | tr '\0' '\n' | sed '/^$/d' | sort -u )
}
# git add must stage exactly what take_list listed before it, or the
# commit would hold something the dry run did not show
add_checked() {
  local want got
  want=$(take_list "$1") || exit 1
  ( cd "$R/frmtmb-wt-$1" && git add -A "${SPEC[@]}" )
  got=$(staged_list "$1") || { echo "wt-$1: git diff-index --cached failed"; exit 1; }
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
  echo "wt-r84-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
rec() { echo "Lane record: dev/$1-findings.md; review: dev/reviews/2026-11-17-$1.md."; }
mrec() { echo "Lane record: dev/$1-findings.md (a measurement lane: dev/ files only, no review)."; }
commit_lane sparsese "frmtmb: REML standard errors above the joint limit by sparse solves" \
"BREAKING: a REML fit with more joint parameters than
options(frmtmb.joint_root_max = ) reads its coefficients' and its
predictions' standard errors by sparse solves in the coordinates of the
design; frm_lp_basis(), frm_curve() and cumulative fitted() follow.
$(rec sparsese)"
commit_lane gpselimit "frmtmb: an approximate gp() at the white limit of its basis" \
"BREAKING: an approximate gp(k = ) whose length scale ran below the
spacing its basis resolves gets a message where the ridge is proved
level; gp(by = ) sub-GPs lose their sdgp and length-scale standard
errors there; the reason \"gpwhite\". $(rec gpselimit)"
commit_lane slopetrig "frmtmb: autoscale on a slope column far from 0; short correlations" \
"BREAKING: the default autoscale engages on a random slope whose column
is far from 0 against its spread, over the rows a block covers; a us()
correlation left short of its maximum is not called a boundary.
$(rec slopetrig)"
commit_lane alapse "frmtmb: asym_laplace() without random effects: no Hessian SE" \
"BREAKING: asym_laplace() and zero_inflated_asym_laplace() fits without
random effects give no Hessian standard error (frmtmb_ald_no_se, reason
\"kink\"); the fit ends at the exact quantile regression solution, with
code 0 only where the maximum is proved global (else frmtmb_ald_local);
all-censored separation; nine more structures refused under Laplace.
$(rec alapse)"
commit_lane nlband "frmtmb: a nonlinear fit starts inside each hole of the data" \
"BREAKING: a nonlinear fit with a pole also starts inside each hole of
the data and keeps a hole fit only where the data place the pole
(frmtmb_nl_pole_refused). $(rec nlband)"
commit_lane quadesc "frmtmb: an abandoned stationary escape keeps its best point" \
"BREAKING: a restart of the escape from a family's declared stationary
point that the optimizer abandons on an error keeps the lowest point
with a finite gradient (escape_track()). $(rec quadesc)"
commit_lane nbrobust "frmtmb: a stable negative-binomial log density" \
"BREAKING: negbinomial(), negbinomial2, nbinom1 and their zero-inflated
and hurdle forms use nb_lpmf(), which keeps its digits and gradient at a
large shape; RTMB-4 and R-1 in dev/upstream-bugs.md. $(rec nbrobust)"
commit_lane squant "frmtmb: quantile_elf(), a smoothed quantile family" \
"New family quantile_elf(): the extended log-F density of qgam, with
curvature in mu, Hessian standard errors and a Laplace maximum with a
group-level term; its compat rows. $(rec squant)"
commit_lane drawdir "frmtmb: predict()'s draw check reads each direction" \
"BREAKING: predict()'s draw check reads each parameter's direction of
the normal law; a fire along variance components only warns and does
not reweight. $(rec drawdir)"
commit_lane ordstarts "frmtmb: ordinal mixtures run no jittered starts by default" \
"BREAKING: no mixture() of ordinal families runs jittered starts by
default; one swapped-component refit for distinct families; the ranking
of the jittered starts reads the degenerate and SE checks.
Lane record: dev/ordstarts-findings.md; reviews:
dev/reviews/2026-11-17-ordstarts.md, dev/reviews/2026-11-18-ordstarts-build.md."
commit_lane cdfname "frmtmb: the CDF slot of frmtmb_family() is cdf" \
"BREAKING: frmtmb_family(lcdf = ) is now cdf =, and the old name is
refused; frmtmb.eam and frmtmb.survival fill the renamed slot; the
rename script dev/cdfname-rename.R. $(rec cdfname)"
commit_lane samplers "dev: brms with stanli and SparseNUTS, measured" \
"The samplers measurement lane: stanli and SparseNUTS against frm_sample()
and brms; upstream entries brms-53, stanli-1, StanEstimators-1 and
SparseNUTS-1. $(mrec samplers)"
commit_lane codesize "dev: the code size of frmtmb, measured" \
"The codesize measurement lane: core R/ size and growth, dead code,
duplication and move-out candidates. $(mrec codesize)"
commit_lane cran "dev: R CMD check --as-cran without NOT_CRAN, measured" \
"The cran measurement lane: the CRAN-tier check time and a tier
proposal. $(mrec cran)"

REL="$R/frmtmb-wt-rel084"
add_checked rel084
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.84.0: eleven lanes and three measures of the round of 2026-11-16

Merges the eleven reviewed lanes (sparsese, gpselimit, slopetrig,
alapse, nlband, quadesc, nbrobust, squant, drawdir, ordstarts, cdfname)
and the dev/ records of three measurement lanes (samplers, codesize,
cran), all on b255350e. frmtmb 0.84.0 (BREAKING). frmtmb.eam 0.14.0 and
frmtmb.survival 0.4.0 need frmtmb 0.84.0 (the cdf slot);
frmtmb.coupling 0.6.2 (compat text).

asym_laplace()'s refusal names quantile_elf() first (merge step 5).

Consolidation: dev/rel084-merge.sh, dev/rel084-news.R,
dev/rel084-release-fix.sh, dev/cdfname-rename.R.

Record: dev/round-20261116.md.

$trailer
EOF
)
ps=(-p "$RELBASE")
for l in $LANES; do ps+=(-p "wt-r84-$l"); done
rel=$(git -C "$REL" commit-tree "$tree" "${ps[@]}" -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel084 now: $(git -C "$REL" log --oneline -1) | left: $(bash "$REL/dev/rel084-list.sh" "$REL" | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1)"

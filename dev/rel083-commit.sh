#!/usr/bin/env bash
# Commit each lane of the round of 2026-11-13 on its branch, then the
# tested release tree as one merge commit with 8159568e and the fourteen
# lanes as parents, and fast-forward main. Modeled on dev/rel082-commit.sh.
# NOT RUN by the consolidation: run it only on the user's OK. It pushes
# nothing.
#
#   bash dev/rel083-commit.sh            # commit
#   bash dev/rel083-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write, no lock taken, and
#                                        # each tree's index sha256 stays
#
# The lanes are on ec8178d3 (frmtmb 0.82.0); the release tree and main
# are on 8159568e (0.82.1, pushed), one commit past it (lane ci082, the
# ubuntu CI fix). Each lane commits its whole tree but the machine output
# that dev/.gitignore lets through: a reviewer's output folder
# (dev/review-*-out*/), a lane's copy of the equivalence survey's output
# (dev/*-equiv-out*/), a reviewer's trial-merge folder
# (dev/review-*-merge/), lane compat's race folders of 0.81.0
# (dev/compat-m-race*/), *.progress and *.done, and this round's two in
# lane hygiene: the Python byte-code cache a reviewer's script left in
# .github/scripts/__pycache__/ and the reviewer's copy of the base
# codemeta.json (dev/reviews/hygiene-review-cm-base.json). Every `git add`
# here excludes those, and the script stops if one would still be taken.
# The release tree carries copies of the lanes' records, reviews and
# evidence scripts (dev/rel083-merge.sh).
#
# main's checkout must be clean: the script stops before any write if it
# is not.
#
# Resuming after a partial run. The lane commits are made one at a time.
# If the script stops after some of them, those lanes' branches are one
# commit past ec8178d3, and the HEAD checks stop a rerun. Do not reset
# them. For each lane already committed, check `git -C <tree> log
# --oneline -2` (its commit, then ec8178d3) and that bash
# dev/rel083-list.sh <tree> lists nothing but the excluded paths. Then
# run the remaining commit_lane calls and the release block by hand, in
# this file's order: the release commit's parents are 8159568e and the
# fourteen wt-r83-<lane> branches, whichever run made them. If the
# release commit exists and main did not move, run only
# `git -C <main checkout> merge --ff-only <release commit>`.
#
# The dry run takes no lock and writes nothing in git: its git calls are
# rev-parse, branch --show-current, ls-files, hash-object (no -w),
# diff-tree between two commits, show, merge-file -p on scratch copies and
# merge-base --is-ancestor, under GIT_OPTIONAL_LOCKS=0; never git status
# or git diff, which refresh an index's stat data. It reads the index's
# sha256 of every tree before and after, and stops if one moved. It
# writes its lists to dev/rel083-log/commit-dry-<tree>.txt.
set -eu
# a failing git or filter inside a pipeline fails the pipeline, so a
# list built from one is never silently short
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=ec8178d3
RELBASE=8159568e
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
LANES="selimit2 quadseam2 slopescale nlstop2 priors2 famend silent2 arspat2 refits parity2 wcheck2 gpfield predraw2 hygiene"
SPEC=(-- . ':(exclude,glob)dev/review-*-out*/**' ':(exclude,glob)dev/*-equiv-out*/**'
      ':(exclude,glob)dev/review-*-merge/**' ':(exclude,glob)dev/compat-m-race*/**'
      ':(exclude,glob)dev/*.progress' ':(exclude,glob)dev/*.done'
      ':(exclude,glob).github/scripts/__pycache__/**'
      ':(exclude)dev/reviews/hygiene-review-cm-base.json')
OUTDIR='^(dev/(review-[^/]*-out[^/]*/|[^/]*-equiv-out[^/]*/|review-[^/]*-merge/|compat-m-race[^/]*/|reviews/hygiene-review-cm-base\.json$)|\.github/scripts/__pycache__/)'
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
blob_moved() {  # paths on stdin; prints those whose content left the index
  local p
  p=$(cat)
  [ -n "$p" ] || return 0
  paste -d '\t' <(printf '%s\n' "$p") <(printf '%s\n' "$p" | git hash-object --stdin-paths) |
    awk -F'\t' 'NR == FNR { split($1, m, " "); blob[$2] = m[2]; next }
                !($1 in blob) || blob[$1] != $2 { print $1 }' \
      <(git ls-files -s -z | tr '\0' '\n') -
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
  for t in frmtmb $(for l in $LANES rel083; do printf 'frmtmb-wt-%s ' "$l"; done); do
    sha256sum "$(git -C "$R/$t" rev-parse --path-format=absolute --git-path index)" ||
      { echo "$t: no index to read" >&2; return 1; }
  done
}
# read before the first git call on any tree, so that the dry run's
# "unchanged" covers every read the script makes
before=$(index_sums) || exit 1
for l in $LANES; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-r83-$l is not on $BASE"; exit 1; }
  [ "$(git -C "$R/frmtmb-wt-$l" branch --show-current)" = "wt-r83-$l" ] ||
    { echo "frmtmb-wt-$l is not on branch wt-r83-$l"; exit 1; }
done
[ "$(git -C "$R/frmtmb-wt-rel083" rev-parse --short=8 HEAD)" = "$RELBASE" ] ||
  { echo "wt-rel083 is not on $RELBASE"; exit 1; }
[ "$(git -C "$R/frmtmb-wt-rel083" branch --show-current)" = wt-rel083 ] ||
  { echo "frmtmb-wt-rel083 is not on branch wt-rel083"; exit 1; }
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$RELBASE" ] ||
  { echo "main is not on $RELBASE"; exit 1; }
[ "$(git -C "$R/frmtmb" branch --show-current)" = main ] ||
  { echo "the main checkout is not on branch main"; exit 1; }
# main must be clean, read as take_list reads a lane (no git status)
mst=$(cd "$R/frmtmb" && {
        git ls-files -z --others --exclude-standard | tr '\0' '\n'
        git ls-files -z --deleted | tr '\0' '\n'
        git ls-files -z --modified | tr '\0' '\n' | sed '/^$/d' | blob_moved
      } | sed '/^$/d' | sort -u) || { echo "main: git ls-files failed"; exit 1; }
[ -z "$mst" ] || { echo "main has uncommitted changes; not fast-forwarding:"; echo "$mst"; exit 1; }
for l in $LANES; do
  [ -f "$R/frmtmb-wt-rel083/dev/reviews/2026-11-13-$l.md" ] ||
    { echo "no review of $l in wt-rel083"; exit 1; }
  [ -f "$R/frmtmb-wt-rel083/dev/$l-findings.md" ] ||
    { echo "no findings of $l in wt-rel083"; exit 1; }
done
for f in dev/round-20261113.md dev/rel083-merge.sh dev/rel083-list.sh \
         dev/rel083-release-fix.sh dev/rel083-news.R dev/rel083-nl-starts.R \
         dev/rel083-fix-at0.R dev/rel083-fix-families.R dev/rel083-fix-fit.R \
         dev/rel083-fix-quadmap.R tests/testthat/test-nl-starts-control.R \
         dev/rel083-log.md; do
  [ -f "$R/frmtmb-wt-rel083/$f" ] || { echo "no $f in wt-rel083"; exit 1; }
done
# no machine output in any commit
for l in $LANES rel083; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | { grep -E "$OUTDIR" || true; })
  [ -z "$bad" ] || { echo "wt-$l: machine output would be committed:"; echo "$bad"; exit 1; }
done
# The drift guard (0.82.0's release review, major 2). A lane that moved
# after the merge would be committed on its branch and made a parent of a
# release commit whose tree lacks the move; git would then count the
# change as merged, and no later merge would bring it in. So each lane
# must take the number of paths the merge took (dev/rel083-merge.sh, its
# taken list), and every path that one lane alone takes must equal the
# release tree's copy, line endings aside, but the paths the consolidation
# edited after the merge (CONS_EDITS, listed in dev/round-20261113.md).
# A single-lane path that 8159568e also changed (lane ci082) cannot equal
# the lane's copy; it is checked as a multi-lane path is.
declare -A MERGED=([selimit2]=66 [quadseam2]=62 [slopescale]=34 [nlstop2]=51
  [priors2]=103 [famend]=82 [silent2]=57 [arspat2]=66 [refits]=85
  [parity2]=85 [wcheck2]=43 [gpfield]=84 [predraw2]=52 [hygiene]=72)
# extensions/frmtmb.spatial/NEWS.md: rebuilt (dev/rel083-news.R);
# R/fit-restart.R: merge step 5; codemeta.json: regenerated at 0.83.0;
# man/par_template.Rd: roxygen output of the merged tree (merge step 5);
# dev/release/run-session.R: pointed at rellib-r22 (merge step 14)
CONS_EDITS=" extensions/frmtmb.spatial/NEWS.md R/fit-restart.R codemeta.json
  man/par_template.Rd dev/release/run-session.R "
# lane:path pairs whose lane lines the consolidation rewrote: the at0
# block (merge step 2), mixture_ord_degenerate_check() (merge step 1),
# the end of frm() (refits x gpfield), the quadrature restart branch
# (merge step 3, next to quadseam2's lines)
MULTI_EXEMPT=" priors2:tests/testthat/test-prior-table-brms.R
  famend:tests/testthat/test-prior-table-brms.R silent2:R/families.R
  famend:R/families.R refits:R/fit.R gpfield:R/fit.R quadseam2:R/fit.R "
RELMOVED=$(git -C "$R/frmtmb-wt-rel083" diff-tree -r --name-only "$BASE" "$RELBASE") ||
  { echo "git diff-tree failed"; exit 1; }
merges_clean() {  # lane path tmpdir: the lane's copy merges into the release copy unchanged
  local rc=0
  git -C "$R/frmtmb-wt-rel083" show "$BASE:$2" 2>/dev/null | tr -d '\r' > "$3/base" || : > "$3/base"
  tr -d '\r' < "$R/frmtmb-wt-rel083/$2" > "$3/rel"
  tr -d '\r' < "$R/frmtmb-wt-$1/$2" > "$3/lane"
  git merge-file -p "$3/rel" "$3/base" "$3/lane" > "$3/out" 2>/dev/null || rc=$?
  [ "$rc" = 0 ] && cmp -s "$3/out" "$3/rel"
}
# When a lane's copy does not merge into the release copy by itself (a
# hunk of another lane sits next to its own, which git merge-file calls
# a conflict), replay the merge as dev/rel083-merge.sh made it: from
# 8159568e's blob, every lane that takes the path, in merge order. The
# release copy must be what the replay gives.
replays_clean() {  # path lanes tmpdir
  local q rc
  git -C "$R/frmtmb-wt-rel083" show "$BASE:$1" 2>/dev/null | tr -d '\r' > "$3/base" || : > "$3/base"
  git -C "$R/frmtmb-wt-rel083" show "$RELBASE:$1" 2>/dev/null | tr -d '\r' > "$3/cur" || : > "$3/cur"
  for q in $LANES; do
    case " $2 " in *" $q "*) ;; *) continue ;; esac
    tr -d '\r' < "$R/frmtmb-wt-$q/$1" > "$3/lane"
    rc=0; git merge-file -q "$3/cur" "$3/base" "$3/lane" 2>/dev/null || rc=$?
    [ "$rc" = 0 ] || return 1
  done
  cmp -s "$3/cur" <(tr -d '\r' < "$R/frmtmb-wt-rel083/$1")
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
    if printf '%s\n' "$RELMOVED" | grep -qxF -- "$p"; then
      merges_clean "$l" "$p" "$tmp" ||
        { echo "DRIFT wt-$l $p: its copy does not merge into the release copy unchanged" >&2; bad=$((bad + 1)); }
    elif [ ! -e "$R/frmtmb-wt-rel083/$p" ] ||
         ! cmp -s <(tr -d '\r' < "$R/frmtmb-wt-$l/$p") <(tr -d '\r' < "$R/frmtmb-wt-rel083/$p"); then
      echo "DRIFT wt-$l $p: differs from the release tree" >&2
      bad=$((bad + 1))
    fi
  done <<< "$single"
  [ "$bad" = 0 ] || { rm -rf "$tmp"; echo "$bad single-lane paths drifted; re-merge them" >&2; return 1; }
  # A path several lanes take: each lane's copy, merged three-way into the
  # release copy over ec8178d3's blob, must give the release copy back, or
  # the lane holds a change the release tree lacks. NEWS files, the
  # backlog and dev/upstream-bugs.md were merged with --union and rebuilt,
  # and the consolidation rewrote lanes' lines in the pairs of MULTI_EXEMPT.
  local multi nmulti=0 nexempt=0 nreplay=0 takers
  multi=$(printf '%s' "$all" |
            awk 'NF { c[$2]++; w[$2] = w[$2] " " $1 } END { for (p in c) if (c[p] > 1) { n = split(w[p], l, " "); for (i = 1; i <= n; i++) print l[i], p } }' |
            sort)
  while read -r l p; do
    [ -n "$p" ] || continue
    case "$p" in NEWS.md|extensions/*/NEWS.md|dev/test-backlog.md|dev/upstream-bugs.md) nexempt=$((nexempt + 1)); continue ;; esac
    case "$MULTI_EXEMPT" in *[[:space:]]"$l:$p"[[:space:]]*) nexempt=$((nexempt + 1)); continue ;; esac
    nmulti=$((nmulti + 1))
    merges_clean "$l" "$p" "$tmp" || {
      takers=$(printf '%s\n' "$multi" | awk -v p="$p" '$2 == p { print $1 }' | tr '\n' ' ')
      if replays_clean "$p" "$takers" "$tmp"; then nreplay=$((nreplay + 1))
      else echo "DRIFT wt-$l $p: its copy does not merge into the release copy unchanged" >&2; bad=$((bad + 1)); fi; }
  done <<< "$multi"
  rm -rf "$tmp"
  [ "$bad" = 0 ] || { echo "$bad multi-lane pairs drifted; re-merge them" >&2; return 1; }
  echo "drift guard: 14 lane counts as merged; $nsingle single-lane paths equal the release tree or merge into it unchanged (consolidation edits aside); $nmulti multi-lane lane-path pairs merge into the release copy unchanged, $nreplay of them by a replay of the merge ($nexempt exempt)"
}
drift_check || exit 1
if $DRY; then
  mkdir -p "$R/frmtmb-wt-rel083/dev/rel083-log"
  for l in $LANES rel083; do
    lst=$(take_list "$l") || exit 1
    printf '%s\n' "$lst" | sed "s#^#   #" \
      > "$R/frmtmb-wt-rel083/dev/rel083-log/commit-dry-$l.txt"
    n=$(grep -c '' "$R/frmtmb-wt-rel083/dev/rel083-log/commit-dry-$l.txt")
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
  echo "wt-r83-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
rec() { echo "Lane record: dev/$1-findings.md; review: dev/reviews/2026-11-13-$1.md."; }
commit_lane selimit2 "frmtmb: the SE check's limit rule per linear predictor" \
"BREAKING: the standard-error check finds coefficients of one linear
predictor that run to a limit together; zero-inflation and hurdle runs;
fitted() NaN standard errors for a row a limit direction reaches; a
gradient along a limit direction; predictions of whitened REML fits
above the joint limit. $(rec selimit2)"
commit_lane quadseam2 "frmtmb: quadrature escapes and the stationary-point escape" \
"BREAKING: quadrature fits on one block per row escape a stop at a
small sd and keep 0.82.0's rule; the escape from a family's declared
stationary point; the latent residual's Laplace warning.
$(rec quadseam2)"
commit_lane slopescale "frmtmb: autoscale of random-effect columns" \
"BREAKING: autoscale rescales a random-effect column whether or not it
is centered; an sd left at its start is said. $(rec slopescale)"
commit_lane nlstop2 "frmtmb: a nonlinear fit beside a pole of its body" \
"BREAKING: a nonlinear fit that stops between the poles of its body
restarts beyond them; the frmtmb_nl_limit warning; ?huber gives the
measured risk of the Laplace approximation. $(rec nlstop2)"
commit_lane priors2 "frmtmb: class theta is a mixture's weights; thetacov" \
"BREAKING: class \"theta\" is brms's mixture weights (dirichlet()), the
internal covariance class is \"thetacov\"; cor with resp refused; bounds
on Intercept, open ends of a prior's support, a bound through a link.
frmtmb.sample reads the theta prior. $(rec priors2)"
commit_lane famend "frmtmb: family ends, starts inside a prior's support" \
"frechet() and gen_extreme_value() at the ends of their parameters; a
cs() effect that closes an ordinal category; a start where a prior's
density is 0 moves inside; BREAKING: OSA residuals of censored counts
and of rows off the CDF. $(rec famend)"
commit_lane silent2 "frmtmb: silent wrong answers of anova, predictions, mixtures" \
"BREAKING: anova() compares the rows each likelihood reads; addition
terms that read no column; group-level designs of predictions; ordinal
mixture thresholds the data do not place; frmtmb.sample's laplace
conditional effects and print order. $(rec silent2)"
commit_lane arspat2 "frmtmb: mi() beside fcor() and sar() by the exact margin" \
"BREAKING: a missing response of fcor() or a registered term leaves the
fit by the exact margin of the observed rows; the margin slot of
frmtmb_register_autocor(); clamped Wald intervals of bounded correlation
parameters. frmtmb.spatial sar() margin. $(rec arspat2)"
commit_lane refits "frmtmb: response coding, refit(), update() and influence()" \
"BREAKING: character ordinal responses refused as brms does; refit() of
bernoulli reads values first; update() keeps the fit's na.action, se
and verbose; influence() keeps the smooth basis; vcov_cluster() reads
the data. $(rec refits)"
commit_lane parity2 "frmtmb: conditional_effects() and other brms surface" \
"BREAKING: conditional_effects() picks brms's default displays and
sample_new_levels; s() of an expression; lf() replaces a parameter;
REML logLik() df; drop_unused_levels; frmtmb.sample loo readers of a
subset and err[k] draws; ce_snl() and latres_err(). $(rec parity2)"
commit_lane wcheck2 "frmtmb: asym_laplace with random effects refused by proof" \
"BREAKING: the weights check with an observation-level effect; a prior
on a dispersion intercept; asym_laplace() and its zero-inflated form
with a random effect in mu under Laplace refused, with a proof, pointing
to frm_sample(). $(rec wcheck2)"
commit_lane gpfield "frmtmb: gp() kernels, draws and end warnings" \
"gp() takes brms's matern kernels; BREAKING: predict() draws a gp()
term's coefficients with the others; exact and approximate gp() checks
at the end of a fit (frmtmb_gp_laplace, frmtmb_gp_basis). $(rec gpfield)"
commit_lane predraw2 "frmtmb: predict() checks its parameter draws" \
"BREAKING: predict() checks its parameter draws against the fit and
reweights them (PSIS), and a steep side's narrowing no longer jumps.
$(rec predraw2)"
commit_lane hygiene "frmtmb: registrations, compat rows, codemeta, coverage plan" \
"frm_compat() rows for REML with mi() and rate() with cens()/trunc();
a registration replaces an earlier one only when provably the same;
codemeta.json and its test; the coverage plan; ode test mocks.
$(rec hygiene)"

REL="$R/frmtmb-wt-rel083"
add_checked rel083
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.83.0: fourteen lanes of the round of 2026-11-13

Merges the fourteen reviewed lanes (selimit2, quadseam2, slopescale,
nlstop2, priors2, famend, silent2, arspat2, refits, parity2, wcheck2,
gpfield, predraw2, hygiene), all on ec8178d3, onto 8159568e (0.82.1).
frmtmb 0.83.0 (BREAKING). frmtmb.sample 0.30.0 and frmtmb.spatial 0.5.0
need frmtmb 0.83.0; frmtmb.latent 0.7.5 needs it for its tests.

frmtmb_control(nl_starts = ) switches off the sign probe and the pole
restart of a nonlinear fit (tests/testthat/test-nl-starts-control.R).

Consolidation: dev/rel083-merge.sh, dev/rel083-fix-*.R,
dev/rel083-news.R, dev/rel083-release-fix.sh, dev/rel083-nl-starts.R.

Record: dev/round-20261113.md.

$trailer
EOF
)
ps=(-p "$RELBASE")
for l in $LANES; do ps+=(-p "wt-r83-$l"); done
rel=$(git -C "$REL" commit-tree "$tree" "${ps[@]}" -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel083 now: $(git -C "$REL" log --oneline -1) | left: $(bash "$REL/dev/rel083-list.sh" "$REL" | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1)"

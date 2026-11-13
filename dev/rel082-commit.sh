#!/usr/bin/env bash
# Commit each lane of the round of 2026-11-11 on its branch, then the
# tested release tree as one merge commit with e8abed0a and the twelve
# lanes as parents, and fast-forward main. Modeled on dev/rel081-commit.sh.
# NOT RUN by the consolidation: run it only on the user's OK. It pushes
# nothing.
#
#   bash dev/rel082-commit.sh            # commit
#   bash dev/rel082-commit.sh --dry-run  # check the HEADs and list what
#                                        # each commit would take; no
#                                        # git write, no lock taken, and
#                                        # each tree's index sha256 stays
#
# The lanes, the release tree and main are all on e8abed0a (frmtmb 0.81.0,
# pushed). Each lane commits its whole tree but the machine output that
# dev/.gitignore lets through: a reviewer's output folder
# (dev/review-*-out*/), a lane's copy of the equivalence survey's output
# (dev/*-equiv-out*/), a reviewer's trial-merge folder
# (dev/review-*-merge/), lane compat's parallel-race output folders of
# 0.81.0 (dev/compat-m-race*/), and this round's two: lane nlstop's
# reviewer copy of the test helpers and a stale test file
# (dev/reviews/nlstop-review-session/, 26 files) and lane names'
# reviewer diffs and profile (dev/reviews/names-review-*.diff, .prof).
# Every `git add` here excludes those,
# and the script stops if one would still be taken. The release tree
# carries copies of the lanes' records, reviews and evidence scripts
# (dev/rel082-merge.sh).
#
# main's checkout must be clean: the script stops before any write if it
# is not.
#
# Resuming after a partial run. The lane commits are made one at a time.
# If the script stops after some of them, those lanes' branches are one
# commit past e8abed0a, and the HEAD checks stop a rerun. Do not reset
# them. For each lane already committed, check `git -C <tree> log
# --oneline -2` (its commit, then e8abed0a) and `git -C <tree> status
# --porcelain` (nothing but the excluded folders). Then run the
# remaining commit_lane calls and the release block by hand, in this
# file's order: the release commit's parents are e8abed0a and the twelve
# wt-r82-<lane> branches, whichever run made them. If the release
# commit exists and main did not move, run only
# `git -C <main checkout> merge --ff-only <release commit>`.
#
# The dry run takes no lock and writes nothing in git: its git calls are
# rev-parse, branch --show-current, status --porcelain, ls-files and
# merge-base --is-ancestor, under GIT_OPTIONAL_LOCKS=0, and it reads the
# index's sha256 of every tree before and after, and stops if one moved.
# It writes its lists to dev/rel082-log/commit-dry-<tree>.txt.
set -eu
# a failing git or filter inside a pipeline fails the pipeline, so a
# list built from one is never silently short
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
BASE=e8abed0a
DRY=false
[ "${1:-}" = --dry-run ] && DRY=true
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
LANES="boundfam quadseam ppdraw wcheck priortab nlstop gpstop ordstart silent arspat names sampfront"
SPEC=(-- . ':(exclude,glob)dev/review-*-out*/**' ':(exclude,glob)dev/*-equiv-out*/**'
      ':(exclude,glob)dev/review-*-merge/**' ':(exclude,glob)dev/compat-m-race*/**'
      ':(exclude,glob)dev/reviews/nlstop-review-session/**'
      ':(exclude,glob)dev/reviews/names-review-*.diff'
      ':(exclude,glob)dev/reviews/names-review-*.prof'
      ':(exclude,glob)dev/*.progress' ':(exclude,glob)dev/*.done')
OUTDIR='^dev/(review-[^/]*-out[^/]*/|[^/]*-equiv-out[^/]*/|review-[^/]*-merge/|compat-m-race[^/]*/|reviews/nlstop-review-session/|reviews/names-review-[^/]*\.(diff|prof)$)'
# the paths `git add -A "${SPEC[@]}"` would take in a tree, read-only:
# every untracked file that is not ignored and every deleted file, as
# git status lists them, and every tracked file whose content differs
# from the index. git status (under GIT_OPTIONAL_LOCKS=0) and git ls-files
# write no index, but they also list a file whose time alone moved (lane
# boundfam's seven extension census copies), which git add does not
# stage. Porcelain git diff would tell them apart but refreshes the
# index's stat data and writes it, whatever GIT_OPTIONAL_LOCKS says (the
# release review's short check, S1). So the content is compared here:
# git hash-object --stdin-paths (no -w; it applies the file's CRLF
# filter by its path) against the blob git ls-files -s holds. It stops
# the script when a listing fails, when git status and git ls-files
# disagree, or when the list is empty (every tree here has work to
# commit). Called as `lst=$(take_list l) || exit 1`, never inside a
# pipeline, so its refusal reaches the script.
take_list() {
  local st a b c tracked moved
  st=$(cd "$R/frmtmb-wt-$1" &&
         git status --porcelain=v1 -z --untracked-files=all --no-renames "${SPEC[@]}" |
         tr '\0' '\n') ||
    { echo "wt-$1: git status failed" >&2; return 1; }
  b=$(printf '%s\n' "$st" | sed '/^$/d' | cut -c4- | sort -u)
  c=$(cd "$R/frmtmb-wt-$1" &&
        git ls-files -z --modified --deleted --others --exclude-standard \
          "${SPEC[@]}" | tr '\0' '\n' | sort -u) ||
    { echo "wt-$1: git ls-files failed" >&2; return 1; }
  [ "$b" = "$c" ] ||
    { echo "wt-$1: git status and git ls-files disagree" >&2; return 1; }
  a=$(printf '%s\n' "$st" |
        awk 'substr($0, 1, 2) == "??" || substr($0, 2, 1) == "D" { print substr($0, 4) }')
  tracked=$(printf '%s\n' "$st" |
              awk 'length($0) > 3 && substr($0, 1, 2) != "??" && substr($0, 2, 1) != "D" { print substr($0, 4) }')
  if [ -n "$tracked" ]; then
    moved=$(cd "$R/frmtmb-wt-$1" &&
              paste -d '\t' <(printf '%s\n' "$tracked") \
                <(printf '%s\n' "$tracked" | git hash-object --stdin-paths) |
              awk -F'\t' 'NR == FNR { split($1, m, " "); blob[$2] = m[2]; next }
                          !($1 in blob) || blob[$1] != $2 { print $1 }' \
                <(git ls-files -s -z | tr '\0' '\n') -) ||
      { echo "wt-$1: git hash-object or git ls-files -s failed" >&2; return 1; }
    a=$(printf '%s\n%s\n' "$a" "$moved")
  fi
  a=$(printf '%s\n' "$a" | sed '/^$/d' | sort -u)
  [ -n "$a" ] ||
    { echo "wt-$1: the list of paths to commit is empty" >&2; return 1; }
  printf '%s\n' "$a"
}
# the sha256 of each tree's index, to show the dry run moved none
index_sums() {
  for t in frmtmb $(for l in $LANES rel082; do printf 'frmtmb-wt-%s ' "$l"; done); do
    sha256sum "$(git -C "$R/$t" rev-parse --path-format=absolute --git-path index)" ||
      { echo "$t: no index to read" >&2; return 1; }
  done
}
# read before the first git call on any tree, so that the dry run's
# "unchanged" covers every read the script makes
before=$(index_sums) || exit 1
for l in $LANES; do
  [ "$(git -C "$R/frmtmb-wt-$l" rev-parse --short=8 HEAD)" = "$BASE" ] ||
    { echo "wt-r82-$l is not on $BASE"; exit 1; }
  [ "$(git -C "$R/frmtmb-wt-$l" branch --show-current)" = "wt-r82-$l" ] ||
    { echo "frmtmb-wt-$l is not on branch wt-r82-$l"; exit 1; }
done
[ "$(git -C "$R/frmtmb-wt-rel082" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "wt-rel082 is not on $BASE"; exit 1; }
[ "$(git -C "$R/frmtmb-wt-rel082" branch --show-current)" = wt-rel082 ] ||
  { echo "frmtmb-wt-rel082 is not on branch wt-rel082"; exit 1; }
[ "$(git -C "$R/frmtmb" rev-parse --short=8 HEAD)" = "$BASE" ] ||
  { echo "main is not on $BASE"; exit 1; }
[ "$(git -C "$R/frmtmb" branch --show-current)" = main ] ||
  { echo "the main checkout is not on branch main"; exit 1; }
mst=$(git -C "$R/frmtmb" status --porcelain) || { echo "main: git status failed"; exit 1; }
[ -z "$mst" ] || { echo "main has uncommitted changes; not fast-forwarding:"; echo "$mst"; exit 1; }
for l in $LANES; do
  [ -f "$R/frmtmb-wt-rel082/dev/reviews/2026-11-11-$l.md" ] ||
    { echo "no review of $l in wt-rel082"; exit 1; }
  [ -f "$R/frmtmb-wt-rel082/dev/$l-findings.md" ] ||
    { echo "no findings of $l in wt-rel082"; exit 1; }
done
for f in dev/round-20261111.md dev/rel082-merge.sh dev/rel082-release-fix.sh \
         dev/rel082-news.R dev/rel082-backlog.R dev/rel082-remerge-one.sh \
         tests/testthat/test-cold-start.R dev/rel082-log.md \
         dev/rel082-remerge-boundfam.sh tests/testthat/test-gp-inner-pd.R \
         tests/testthat/test-quad-escape-record.R; do
  [ -f "$R/frmtmb-wt-rel082/$f" ] || { echo "no $f in wt-rel082"; exit 1; }
done
# no machine output folder in any commit
for l in $LANES rel082; do
  lst=$(take_list "$l") || exit 1
  bad=$(printf '%s\n' "$lst" | { grep -E "$OUTDIR" || true; })
  [ -z "$bad" ] || { echo "wt-$l: machine output would be committed:"; echo "$bad"; exit 1; }
done
# The drift guard (release review, major 2). A lane that moved after the
# merge would be committed on its branch and made a parent of a release
# commit whose tree lacks the move; git would then count the change as
# merged, and no later merge would bring it in. So each lane must take the
# number of paths the merge took (dev/rel082-merge.sh, its taken list),
# and every path that one lane alone takes must equal the release tree's
# copy, line endings aside, but the paths the consolidation edited after
# the merge, which are listed here and in dev/round-20261111.md. boundfam:
# 50 paths at the merge, 58 after the re-merge of its consolidation fixes,
# 61 after the re-merge of its OpenBLAS fix (dev/rel082-remerge-boundfam.sh,
# run twice; 7 census copies whose time alone moved are not taken).
declare -A MERGED=([boundfam]=61 [quadseam]=59 [ppdraw]=49 [wcheck]=33
  [priortab]=40 [nlstop]=58 [gpstop]=48 [ordstart]=44 [silent]=73
  [arspat]=66 [names]=78 [sampfront]=62)
CONS_EDITS=" R/autoscale.R R/gp-start.R tests/testthat/test-predict-draw-edge.R
  extensions/frmtmb.latent/NEWS.md extensions/frmtmb.spatial/NEWS.md
  extensions/frmtmb.survival/NEWS.md "
# lane:path pairs whose lane lines the consolidation rewrote: the three
# cold-start variables of R/fit.R (merge step 2) and ppdraw's get0()
# calls in R/predict-brms.R (the release review's minor 5)
MULTI_EXEMPT=" nlstop:R/fit.R gpstop:R/fit.R ordstart:R/fit.R
  ppdraw:R/predict-brms.R "
drift_check() {
  local l lst n all="" single p bad=0 nsingle=0
  for l in $LANES; do
    lst=$(take_list "$l") || return 1
    n=$(printf '%s\n' "$lst" | grep -c .)
    [ "$n" = "${MERGED[$l]:-none}" ] ||
      { echo "wt-$l: takes $n paths, the merge took ${MERGED[$l]:-none}" >&2; return 1; }
    all="$all$(printf '%s\n' "$lst" | sed "s#^#$l #")"$'\n'
  done
  single=$(printf '%s' "$all" |
             awk 'NF { c[$2]++; w[$2] = $1 } END { for (p in c) if (c[p] == 1) print w[p], p }' |
             sort)
  while read -r l p; do
    [ -n "$p" ] || continue
    nsingle=$((nsingle + 1))
    case "$CONS_EDITS" in *[[:space:]]"$p"[[:space:]]*) continue ;; esac
    if [ ! -e "$R/frmtmb-wt-rel082/$p" ] ||
         ! cmp -s <(tr -d '\r' < "$R/frmtmb-wt-$l/$p") <(tr -d '\r' < "$R/frmtmb-wt-rel082/$p"); then
      echo "DRIFT wt-$l $p: differs from the release tree" >&2
      bad=$((bad + 1))
    fi
  done <<< "$single"
  [ "$bad" = 0 ] || { echo "$bad single-lane paths drifted; re-merge them" >&2; return 1; }
  # A path several lanes take: each lane's copy, merged three-way into the
  # release copy over e8abed0a's blob, must give the release copy back, or
  # the lane holds a change the release tree lacks (the release review's
  # short check). Read-only: git show and git merge-file -p into a scratch
  # folder. NEWS files and the backlog are rebuilt, not merged, and the
  # consolidation rewrote lanes' lines in two code files (MULTI_EXEMPT).
  local multi tmp nmulti=0 nexempt=0 rc
  multi=$(printf '%s' "$all" |
            awk 'NF { c[$2]++; w[$2] = w[$2] " " $1 } END { for (p in c) if (c[p] > 1) { n = split(w[p], l, " "); for (i = 1; i <= n; i++) print l[i], p } }' |
            sort)
  tmp=$(mktemp -d) || return 1
  while read -r l p; do
    [ -n "$p" ] || continue
    case "$p" in NEWS.md|extensions/*/NEWS.md|dev/test-backlog.md) nexempt=$((nexempt + 1)); continue ;; esac
    case "$MULTI_EXEMPT" in *[[:space:]]"$l:$p"[[:space:]]*) nexempt=$((nexempt + 1)); continue ;; esac
    nmulti=$((nmulti + 1))
    git -C "$R/frmtmb-wt-rel082" show "$BASE:$p" > "$tmp/base" 2>/dev/null || : > "$tmp/base"
    tr -d '\r' < "$R/frmtmb-wt-rel082/$p" > "$tmp/rel"
    tr -d '\r' < "$R/frmtmb-wt-$l/$p" > "$tmp/lane"
    rc=0
    git merge-file -p "$tmp/rel" "$tmp/base" "$tmp/lane" > "$tmp/out" 2>/dev/null || rc=$?
    if [ "$rc" != 0 ] || ! cmp -s "$tmp/out" "$tmp/rel"; then
      echo "DRIFT wt-$l $p: its copy does not merge into the release copy unchanged" >&2
      bad=$((bad + 1))
    fi
  done <<< "$multi"
  rm -rf "$tmp"
  [ "$bad" = 0 ] || { echo "$bad multi-lane pairs drifted; re-merge them" >&2; return 1; }
  echo "drift guard: 12 lane counts as merged; $nsingle single-lane paths equal the release tree (consolidation edits aside); $nmulti multi-lane lane-path pairs merge into the release copy unchanged ($nexempt exempt)"
}
drift_check || exit 1
if $DRY; then
  mkdir -p "$R/frmtmb-wt-rel082/dev/rel082-log"
  for l in $LANES rel082; do
    lst=$(take_list "$l") || exit 1
    printf '%s\n' "$lst" | sed "s#^#   #" \
      > "$R/frmtmb-wt-rel082/dev/rel082-log/commit-dry-$l.txt"
    n=$(grep -c '' "$R/frmtmb-wt-rel082/dev/rel082-log/commit-dry-$l.txt")
    st=$(cd "$R/frmtmb-wt-$l" && git status --porcelain --untracked-files=normal) ||
      { echo "wt-$l: git status failed"; exit 1; }
    left=$(printf '%s\n' "$st" | { grep -E '^\?\? dev/((review-[^/]*-out[^/]*|[^/]*-equiv-out[^/]*|review-[^/]*-merge|compat-m-race[^/]*|reviews/nlstop-review-session)/|reviews/names-review-[^/]*\.(diff|prof))$' || true; } |
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
  echo "wt-r82-$n: $(git -C "$R/frmtmb-wt-$n" log --oneline -1) | $(git -C "$R/frmtmb-wt-$n" show --stat --format= HEAD | tail -1)"
}
rec() { echo "Lane record: dev/$1-findings.md; review: dev/reviews/2026-11-11-$1.md."; }
commit_lane boundfam "frmtmb: a parameter at the end of its range, told by the SE check" \
"BREAKING: the standard-error check tells a coefficient run to the end
of its range from a flat ridge: a family parameter at its limit gets a
message, a mean coefficient a warning, a separating probability slope
warns; fitted() keeps its standard errors at a held limit and where a
limit direction barely moves it; the quiet limit message has the
subclass frmtmb_limit_fit, which the nonlinear-start restart counts.
$(rec boundfam)"
commit_lane quadseam "frmtmb: quadrature fits settle their anchor; headroom fallback" \
"BREAKING: a quadrature fit re-anchors its rule until the optimum stays,
and warns when it does not settle; the escape order under quadrature;
a finite-difference fallback for the gradient and convergence headroom.
$(rec quadseam)"
commit_lane ppdraw "frmtmb: predict() draws at an edge from a per-side law" \
"BREAKING: predict() draws a parameter at the edge of its range from a
law on each side, with the conditional sd as reference, and reports
the rows its simulation leaves NA. $(rec ppdraw)"
commit_lane wcheck "frmtmb: weighted likelihoods with no maximum, refused by proof" \
"BREAKING: the weights check refuses a model only where it proves the
objective unbounded, with the Jacobian of the shape's link; a model it
cannot prove unbounded fits as on 0.81.0. $(rec wcheck)"
commit_lane priortab "frmtmb: prior tables with bounds, a support box, brms's rows" \
"BREAKING: prior tables carry lb and ub; a table passed back stays ML;
cor rows per group; cor with dpar or nlpar and rescor with coef or group
refused; resp by its Stan name; a support box only where the density is
finite at the end. frmtmb.sample's prior messages. $(rec priortab)"
commit_lane nlstop "frmtmb: a restart at the prior's location; nl limits said" \
"BREAKING: on the user's cold start, a fit with a free dispersion whose
coefficient lies far from its prior's location restarts there and keeps
the higher mode (frmtmb_prior_multistart); a nonlinear body on a ridge
to a limit says so. huber() with group terms documented. $(rec nlstop)"
commit_lane gpstop "frmtmb: a gp() in a nonlinear body is read a second time" \
"BREAKING: an exact gp() in a parameter that a nonlinear body multiplies
is linearized again for its field on the user's cold start, and a
white-noise field is refused as a maximum. $(rec gpstop)"
commit_lane ordstart "frmtmb: jittered starts for ordinal mixtures" \
"BREAKING: a mixture of ordinal families without random effects runs 10
jittered starts by default (frmtmb_control(starts = )), and a step
function yields to the best converged proper start, with the warning
frmtmb_mixture_step_start. $(rec ordstart)"
commit_lane silent "frmtmb: silent wrong answers of DHARMa, nobs, drop1, cox, hmm" \
"BREAKING: dharma_residuals() of cens() rows; nobs(), df.residual() and
anova() of mv and mi() fits; drop1() and add1() of special terms;
frmtmb.survival cox() upper knot (BREAKING); frmtmb.latent hmm() keeps
its state family's stationary point. $(rec silent)"
commit_lane arspat "frmtmb: se() and mi() beside residual correlation terms" \
"ar(), ma() and arma() with se(); mi() beside fcor() and a registered
term such as sar(); a cluster sandwich; interior time gaps reported;
frmtmb.spatial sar() with complex eigenvalues. $(rec arspat)"
commit_lane names "frmtmb: brms names and layouts of variables, ranef(), coef()" \
"BREAKING: variables() and every draws accessor in brms's order;
ranef() and coef() one entry per grouping factor; brms's names for
smooth and gp() weights, condar states and K categories. $(rec names)"
commit_lane sampfront "frmtmb.sample: subsets of the fitted rows; laplace draws at 0" \
"BREAKING: the draws readers take newdata whose rows are fitted rows;
laplace draws of a prediction that reads no group-level value;
model_weights(), pp_average(), prior_samples(). Core exports
integrated_design_zero(). $(rec sampfront)"

REL="$R/frmtmb-wt-rel082"
add_checked rel082
tree=$(git -C "$REL" write-tree)
msg=$(cat <<EOF
Release frmtmb 0.82.0: twelve lanes of the round of 2026-11-11

Merges the twelve reviewed lanes (boundfam, quadseam, ppdraw, wcheck,
priortab, nlstop, gpstop, ordstart, silent, arspat, names, sampfront),
all on e8abed0a, onto e8abed0a. frmtmb 0.82.0 (BREAKING). frmtmb.sample
0.29.0 and frmtmb.spatial 0.4.0 need frmtmb 0.82.0; frmtmb.survival
0.3.0 (BREAKING knots); frmtmb.latent 0.7.4; frmtmb.eam 0.13.2 and
frmtmb.ode 0.7.4 need frmtmb 0.82.0 for their tests.

One cold-start rule for the restart at a prior's location, the second
reading of a gp() body and the mixture starts; autoscale_choose()
replays the messages of the fit it keeps (tests/testthat/test-cold-start.R).

Consolidation: dev/rel082-merge.sh, dev/rel082-news.R,
dev/rel082-backlog.R, dev/rel082-release-fix.sh; after the release review
(dev/reviews/2026-11-12-rel082.md): dev/rel082-remerge-boundfam.sh (lane
boundfam's consolidation fixes), the restart record kept through later
escapes, gp_chol_pd(), the drift guard in this commit's script.

Record: dev/round-20261111.md.

$trailer
EOF
)
ps=(-p "$BASE")
for l in $LANES; do ps+=(-p "wt-r82-$l"); done
rel=$(git -C "$REL" commit-tree "$tree" "${ps[@]}" -m "$msg")
echo "release commit: $rel"
git -C "$REL" reset -q "$rel"
echo "wt-rel082 now: $(git -C "$REL" log --oneline -1) | left: $(git -C "$REL" status --porcelain | wc -l)"
git -C "$R/frmtmb" merge --ff-only -q "$rel"
echo "main now: $(git -C "$R/frmtmb" log --oneline -1) | status: $(git -C "$R/frmtmb" status --porcelain | wc -l)"

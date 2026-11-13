#!/usr/bin/env bash
# Commit lane ci082 as frmtmb 0.82.1 on branch wt-r83-ci082 and
# fast-forward main to it. A patch on ec8178d3 (0.82.0) for the red ubuntu
# CI of 0.82.0. Modeled on dev/rel082-commit.sh. It pushes nothing.
# NOT RUN without the user's OK.
#
#   bash dev/rel0821-commit.sh --manifest  # after the last short check:
#                                          # record each path the commit
#                                          # takes and its blob hash
#   bash dev/rel0821-commit.sh --dry-run   # check HEADs, list the paths,
#                                          # compare with the manifest; no
#                                          # git write, no lock, index
#                                          # sha256 checked unchanged
#   bash dev/rel0821-commit.sh             # bump, commit, fast-forward
#
# The drift guard: the commit takes exactly the paths in the manifest,
# each with the content it had when the manifest was written, plus the
# version bump in DESCRIPTION and the NEWS heading. Anything else stops
# the script before any write.
set -eu
set -o pipefail
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
T="$R/frmtmb-wt-ci082"
M="$R/frmtmb"
BASE=ec8178d3
BRANCH=wt-r83-ci082
MAN="$T/dev/rel0821-manifest.txt"
MODE=commit
case "${1:-}" in
  --dry-run) MODE=dry ;;
  --manifest) MODE=manifest ;;
  "") ;;
  *) echo "unknown argument: $1"; exit 1 ;;
esac
trailer="Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"

# Paths `git add -A` would take, read-only (as dev/rel082-commit.sh):
# untracked-not-ignored and deleted files from git status, and tracked
# files whose content differs from the index by git hash-object (git diff
# would rewrite the index's stat data).
take_list() {
  local st a tracked moved
  st=$(cd "$T" && git status --porcelain=v1 -z --untracked-files=all --no-renames |
         tr '\0' '\n') || { echo "git status failed" >&2; return 1; }
  a=$(printf '%s\n' "$st" |
        awk 'substr($0, 1, 2) == "??" || substr($0, 2, 1) == "D" { print substr($0, 4) }')
  tracked=$(printf '%s\n' "$st" |
              awk 'length($0) > 3 && substr($0, 1, 2) != "??" && substr($0, 2, 1) != "D" { print substr($0, 4) }')
  if [ -n "$tracked" ]; then
    moved=$(cd "$T" &&
              paste -d '\t' <(printf '%s\n' "$tracked") \
                <(printf '%s\n' "$tracked" | git hash-object --stdin-paths) |
              awk -F'\t' 'NR == FNR { split($1, m, " "); blob[$2] = m[2]; next }
                          !($1 in blob) || blob[$1] != $2 { print $1 }' \
                <(git ls-files -s -z | tr '\0' '\n') -) ||
      { echo "git hash-object or git ls-files -s failed" >&2; return 1; }
    a=$(printf '%s\n%s\n' "$a" "$moved")
  fi
  a=$(printf '%s\n' "$a" | sed '/^$/d' | sort -u)
  [ -n "$a" ] || { echo "the list of paths to commit is empty" >&2; return 1; }
  printf '%s\n' "$a"
}
# path<TAB>blob for each listed path (deleted paths get "deleted")
hash_list() {
  local p
  while IFS= read -r p; do
    if [ -e "$T/$p" ]; then
      printf '%s\t%s\n' "$p" "$(cd "$T" && git hash-object -- "$p")"
    else
      printf '%s\tdeleted\n' "$p"
    fi
  done
}
# drop this script and its manifest from a path<TAB>blob list: both are
# written by this release step and so cannot be hashed into the manifest
own_out() {
  awk -F'\t' '$1 != "dev/rel0821-manifest.txt" && $1 != "dev/rel0821-commit.sh"'
}
index_sums() {
  for t in "$M" "$T"; do
    sha256sum "$(git -C "$t" rev-parse --path-format=absolute --git-path index)" ||
      { echo "$t: no index to read" >&2; return 1; }
  done
}

before=$(index_sums) || exit 1
[ "$(git -C "$T" rev-parse --short=8 HEAD)" = "$BASE" ] || { echo "$BRANCH is not on $BASE"; exit 1; }
[ "$(git -C "$T" branch --show-current)" = "$BRANCH" ] || { echo "frmtmb-wt-ci082 is not on $BRANCH"; exit 1; }
[ "$(git -C "$M" rev-parse --short=8 HEAD)" = "$BASE" ] || { echo "main is not on $BASE"; exit 1; }
[ "$(git -C "$M" branch --show-current)" = main ] || { echo "the main checkout is not on main"; exit 1; }
mst=$(git -C "$M" status --porcelain) || { echo "main: git status failed"; exit 1; }
[ -z "$mst" ] || { echo "main has uncommitted changes:"; echo "$mst"; exit 1; }
grep -qE '^Version: 0\.82\.0\r?$' "$T/DESCRIPTION" || { echo "DESCRIPTION is not at 0.82.0"; exit 1; }
[ "$(head -1 "$T/NEWS.md" | tr -d '\r')" = "# frmtmb (development version)" ] ||
  { echo "NEWS.md does not open with the development-version heading"; exit 1; }
[ "$(grep -c '^# frmtmb (development version)' "$T/NEWS.md")" = 1 ] ||
  { echo "NEWS.md has more than one development-version heading"; exit 1; }
for f in dev/ci082-findings.md dev/reviews/2026-11-13-ci082.md dev/reviews/2026-11-13-ci082-short.md; do
  [ -f "$T/$f" ] || { echo "no $f"; exit 1; }
done

lst=$(take_list) || exit 1
cur=$(printf '%s\n' "$lst" | hash_list)

if [ "$MODE" = manifest ]; then
  # the manifest itself and this script are taken too; they are written
  # here and so not hashed into the manifest
  printf '%s\n' "$cur" | own_out > "$MAN"
  echo "manifest: $(grep -c '' "$MAN") paths written to dev/rel0821-manifest.txt"
  after=$(index_sums) || exit 1
  [ "$before" = "$after" ] || { echo "AN INDEX MOVED"; exit 1; }
  exit 0
fi

[ -f "$MAN" ] || { echo "no manifest; run --manifest after the last short check"; exit 1; }
want=$(own_out < "$MAN" | sort)
got=$(printf '%s\n' "$cur" | own_out | sort)
if [ "$want" != "$got" ]; then
  echo "DRIFT: the tree differs from the manifest:"
  diff <(printf '%s\n' "$want") <(printf '%s\n' "$got") || true
  exit 1
fi
echo "drift guard: $(printf '%s\n' "$want" | grep -c .) paths equal the manifest"

if [ "$MODE" = dry ]; then
  printf '%s\n' "$lst" | sed 's#^#   #'
  echo "plus at commit: DESCRIPTION 0.82.0 -> 0.82.1; NEWS heading -> # frmtmb 0.82.1"
  echo "main: fast-forward possible: $(git -C "$M" merge-base --is-ancestor "$BASE" HEAD && echo yes)"
  after=$(index_sums) || exit 1
  [ "$before" = "$after" ] || { echo "AN INDEX MOVED IN THE DRY RUN"; exit 1; }
  echo "index sha256 of main and wt-ci082 unchanged"
  echo "DRY RUN DONE (no git write, no lock)"
  exit 0
fi

# the bump, keeping each file's line endings
perl -pi -e 's/^Version: 0\.82\.0(\r?)$/Version: 0.82.1$1/' "$T/DESCRIPTION"
perl -pi -e 's/^# frmtmb \(development version\)(\r?)$/# frmtmb 0.82.1$1/ if $. == 1' "$T/NEWS.md"
grep -qE '^Version: 0\.82\.1\r?$' "$T/DESCRIPTION" || { echo "version bump failed"; exit 1; }
[ "$(head -1 "$T/NEWS.md" | tr -d '\r')" = "# frmtmb 0.82.1" ] || { echo "NEWS heading failed"; exit 1; }

want_staged=$( (printf '%s\n' "$lst"; echo DESCRIPTION) | sort -u)
( cd "$T" && git add -A )
got_staged=$(cd "$T" && git diff --cached --name-only --no-renames -z | tr '\0' '\n' | sort -u)
[ "$want_staged" = "$got_staged" ] ||
  { echo "the staged paths differ from the list; nothing committed"; diff <(echo "$want_staged") <(echo "$got_staged") || true; exit 1; }

git -C "$T" commit -q -F - <<EOF
Release frmtmb 0.82.1: the red ubuntu CI of 0.82.0

R-CMD-check on 0.82.0 (ec8178d3) failed on ubuntu release and devel
under OpenBLAS. confint_varcorr() and summary() give a covariance
component that its structure fixes, inside the correlation clamp, zero
width on every platform, decided by the component's values at points far
apart in the structure's parameters, not by the size of a rounding
residue; estimated components are as in 0.82.0. The jittered starts of an ordinal mixture restore every
RNGkind() element, including R-devel's binom.kind. Two tests that pinned
one of several valid endings now accept each valid ending by its text.

Record: dev/ci082-findings.md; reviews: dev/reviews/2026-11-13-ci082.md,
dev/reviews/2026-11-13-ci082-short.md.

$trailer
EOF
echo "$BRANCH: $(git -C "$T" log --oneline -1) | $(git -C "$T" show --stat --format= HEAD | tail -1)"
echo "wt-ci082 left: $(git -C "$T" status --porcelain | wc -l)"
git -C "$M" merge --ff-only -q "$BRANCH"
echo "main now: $(git -C "$M" log --oneline -1) | status: $(git -C "$M" status --porcelain | wc -l)"

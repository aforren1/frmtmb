#!/usr/bin/env bash
# Write dev/spatial-out/AB-final.diff: the stage tree (parts A and B on
# the base commit) against the base commit, one git-style section per
# changed, new or removed file, with a/<path> and b/<path> headers so
# that `git apply` on the base commit takes it. Read-only on git.
#
#   bash dev/spatial-ab-diff.sh <stage tree>
set -e
cd "$(dirname "$0")/.."
ab="$1"
out="$PWD/dev/spatial-out/AB-final.diff"
base="$(dirname "$ab")/spatial-base"
rm -rf "$base"; mkdir -p "$base"
git archive 19b162b9 | tar -x -C "$base"
# the files a build leaves in the stage tree are not part of the change
skip='^(dev/|src/.*\.(o|dll|so)$|src/symbols\.rds$)'
{ (cd "$base" && find . -type f); (cd "$ab" && find . -type f); } |
  sed 's|^\./||' | grep -Ev "$skip" | sort -u > "$base.files"
: > "$out"
while read -r f; do
  a="$base/$f"; b="$ab/$f"
  if [ -f "$a" ] && [ -f "$b" ] && cmp -s "$a" "$b"; then continue; fi
  [ -f "$a" ] || a=/dev/null
  [ -f "$b" ] || b=/dev/null
  git -c core.safecrlf=false diff --no-index --binary "$a" "$b" |
    sed -e "1s|^diff --git .*|diff --git a/$f b/$f|" \
        -e "s|^--- a/.*|--- a/$f|" -e "s|^+++ b/.*|+++ b/$f|" >> "$out" ||
    true
done < "$base.files"
echo "$out: $(grep -c '^diff --git' "$out") files"

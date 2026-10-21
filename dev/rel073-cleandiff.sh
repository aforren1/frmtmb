#!/usr/bin/env bash
# Compare the release tree with a merge run from scratch: every tracked
# and untracked (not ignored) file of the tree against the same path in
# <clean>, a copy made with `git archive 4854e425` plus the
# dev/rel073-* inputs, into which `REL073_ROOT=<clean>
# REL073_GIT=<release tree> bash <clean>/dev/rel073-merge.sh` merged
# the lanes. Prints DIFF and ONLY-REL lines (line endings ignored).
#   bash dev/rel073-cleandiff.sh <clean>
cd /c/Users/adf44/source/r/frmtmb-wt-release || exit 1
c="$1"; n=0
while read -r f; do
  [ -f "$f" ] || continue
  n=$((n + 1))
  if [ ! -f "$c/$f" ]; then echo "ONLY-REL $f"
  elif ! diff -q --strip-trailing-cr "$f" "$c/$f" > /dev/null; then
    echo "DIFF $f"
  fi
done < <({ git ls-files; git ls-files --others --exclude-standard; } | sort -u)
echo "compared $n files"

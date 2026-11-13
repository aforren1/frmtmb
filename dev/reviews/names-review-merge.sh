#!/usr/bin/env bash
# Review of lane names: three-way merge of every file the lane shares with
# another lane's worktree, LF-normalized, base e8abed0a. Read-only on the
# worktrees; output under dev/reviews/names-review-merge/.
export GIT_OPTIONAL_LOCKS=0
R=/c/Users/adf44/source/r
N=$R/frmtmb-wt-names
out=$N/dev/reviews/names-review-merge
rm -rf "$out"; mkdir -p "$out"
cd $N || exit 1
mine=$( { git diff --name-only e8abed0a; git ls-files --others --exclude-standard; } | grep -v '^dev/names-' | sort -u)
for w in $R/frmtmb-wt-*; do
  ln=$(basename $w); ln=${ln#frmtmb-wt-}
  [ "$ln" = names ] && continue
  theirs=$(cd $w && { git diff --name-only e8abed0a; git ls-files --others --exclude-standard; } 2>/dev/null | sort -u)
  for f in $(comm -12 <(echo "$mine") <(echo "$theirs")); do
    case "$f" in man/*) continue;; esac
    k=$(echo "$f" | tr '/' '_')
    git show e8abed0a:"$f" 2>/dev/null | tr -d '\r' > "$out/$ln.$k.base"
    tr -d '\r' < "$N/$f" > "$out/$ln.$k.ours"
    tr -d '\r' < "$w/$f" > "$out/$ln.$k.theirs"
    cp "$out/$ln.$k.ours" "$out/$ln.$k.merged"
    git merge-file --diff3 -L names -L base -L "$ln" "$out/$ln.$k.merged" "$out/$ln.$k.base" "$out/$ln.$k.theirs"
    rc=$?
    echo "MERGE $ln $f conflicts=$rc"
  done
done

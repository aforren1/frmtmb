#!/usr/bin/env bash
# Write dev/sarfcor-out/A.diff (the base commit to the tree with part A
# alone) and dev/sarfcor-out/B.diff (that tree to the final one), one
# git-style section per changed, new or removed file, with a/<path> and
# b/<path> headers, so that `git apply` takes A.diff on 18b536fc and
# B.diff after it. Read-only on git (git archive and diff --no-index).
#
# The part-A tree is the base plus the stage-A overlay, the copy of
# every file part A changed, kept in the session scratchpad while part B
# was written (its man/ and NAMESPACE are regenerated here by
# roxygenise). The final tree is the base plus every changed or new file
# of the worktree. dev/ is left out but for the three dev files the
# change edits (dev/core-boundary.md, dev/release/run-session.R,
# dev/upstream-bugs.md). Line endings are normalized to LF.
#
#   bash dev/sarfcor-ab-diff.sh
set -e
WT=C:/Users/adf44/source/r/frmtmb-wt-sarfcor
SP=C:/Users/adf44/AppData/Local/Temp/12/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad
STAGE="$SP/sarfcor-stageA"
T="$SP/sarfcor-trees"
OUT="$WT/dev/sarfcor-out"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
keep_dev='^dev/(core-boundary\.md|release/run-session\.R|upstream-bugs\.md)$'
rm -rf "$T"; mkdir -p "$T/base" "$T/A" "$T/final"
git -C "$WT" archive 18b536fc | tar -x -C "$T/base"
cp -r "$T/base/." "$T/A/"
cp -r "$T/base/." "$T/final/"
cp -r "$STAGE/." "$T/A/"
cd "$WT"
{ git diff --name-only; git ls-files --others --exclude-standard; } |
  sort -u | while read -r f; do
  case "$f" in dev/*) echo "$f" | grep -Eq "$keep_dev" || continue;; esac
  if [ -f "$f" ]; then
    mkdir -p "$T/final/$(dirname "$f")"; cp "$f" "$T/final/$f"
  else
    rm -f "$T/final/$f"
  fi
done
# the A tree's generated files
"$RS" "$WT/dev/sarfcor-rox-tree.R" "$T/A" | grep "ROX DONE"
# a second pass: the first can write NAMESPACE before it reads every
# import tag (frmtmb.spatial lost one import in the first pass)
"$RS" "$WT/dev/sarfcor-rox-tree.R" "$T/A" | grep "ROX DONE"
# normalize text files outside dev/ to LF
for tr in base A final; do
  (cd "$T/$tr" && find . -type f ! -path './dev/*' -print0 |
     xargs -0 grep -Il . | while read -r f; do sed -i 's/\r$//' "$f"; done)
  for f in dev/core-boundary.md dev/release/run-session.R dev/upstream-bugs.md; do
    [ -f "$T/$tr/$f" ] && sed -i 's/\r$//' "$T/$tr/$f"
  done
done
diff_trees() {
  local a="$1" b="$2" out="$3"
  : > "$out"
  { (cd "$a" && find . -type f); (cd "$b" && find . -type f); } |
    sed 's|^\./||' | sort -u |
    grep -Ev '^dev/' > "$T/files.lst" || true
  printf '%s\n' dev/core-boundary.md dev/release/run-session.R \
    dev/upstream-bugs.md >> "$T/files.lst"
  while read -r f; do
    x="$a/$f"; y="$b/$f"
    if [ -f "$x" ] && [ -f "$y" ] && cmp -s "$x" "$y"; then continue; fi
    [ -f "$x" ] || x=/dev/null
    [ -f "$y" ] || y=/dev/null
    git -c core.safecrlf=false diff --no-index --binary "$x" "$y" |
      sed -e "1s|^diff --git .*|diff --git a/$f b/$f|" \
          -e "s|^--- a/.*|--- a/$f|" -e "s|^+++ b/.*|+++ b/$f|" >> "$out" ||
      true
  done < "$T/files.lst"
  echo "$out: $(grep -c '^diff --git' "$out") files"
}
diff_trees "$T/base" "$T/A" "$OUT/A.diff"
diff_trees "$T/A" "$T/final" "$OUT/B.diff"

#!/usr/bin/env bash
# Merge lane nanse onto the 0.68.0 release tree, following the nanse
# review's "Final merge recipe" (dev/reviews/2026-10-06-nanse.md), with
# the lane's CURRENT R/se-check.R and tests/testthat/test-se-check.R
# (its punch round 2b, RB4). Rerunnable: the first run saves the
# release's pre-merge copy of every file the merge touches in
# dev/rel068-log/nanse-pre/, and every run starts from that copy, so a
# replacement file from the lane can be merged again by rerunning.
#
#   bash dev/rel068-nanse-merge.sh
# No git write: git merge-file works on files, not on the index.
set -u
LANE=/c/Users/adf44/source/r/frmtmb-wt-nanse
REL=/c/Users/adf44/source/r/frmtmb-wt-release
PRE=$REL/dev/rel068-log/nanse-pre
W=$REL/dev/rel068-log/nanse-merge-work
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
BASE=9e902909
cd "$REL" || exit 1
# the lane's changed tracked files, without its generated Rd files and
# the two NEWS files, which dev/rel068-news-merge.R writes
files=$(git -C "$LANE" diff "$BASE" --name-only 2>/dev/null |
          grep -vE '^man/.*\.Rd$|^NEWS\.md$|^extensions/frmtmb\.spline/NEWS\.md$')
# the files the resolver edits beyond those
extra="tests/testthat/test-gp-by.R tests/testthat/test-ordinal-mixture.R
extensions/frmtmb.spline/tests/testthat/test-difference.R
extensions/frmtmb.learn/DESCRIPTION extensions/frmtmb.learn/NEWS.md
extensions/frmtmb.sample/R/sample.R"
if [ ! -d "$PRE" ]; then
  mkdir -p "$PRE"
  for f in $files $extra; do
    mkdir -p "$PRE/$(dirname "$f")"; cp -p "$f" "$PRE/$f"
  done
  echo "saved $(find "$PRE" -type f | wc -l) pre-merge files"
fi
for f in $files $extra; do cp -p "$PRE/$f" "$f"; done
rm -rf "$W"; mkdir -p "$W"
n_clean=0; n_conf=0
for f in $files; do
  git -C "$LANE" show "$BASE:$f" | tr -d '\r' > "$W/base"
  tr -d '\r' < "$LANE/$f" > "$W/lane"
  tr -d '\r' < "$f" > "$W/rel"
  if [ "$f" = extensions/frmtmb.ode/tests/testthat/test-ode-nlf.R ]; then
    # step 6: the release's version (lane fixes' flat warning)
    cp "$W/rel" "$f"; echo "   release $f"; continue
  fi
  git merge-file -p -L release -L base -L lane "$W/rel" "$W/base" \
    "$W/lane" > "$W/out"
  rc=$?
  if [ "$rc" -eq 0 ]; then
    n_clean=$((n_clean + 1))
  else
    n_conf=$((n_conf + 1)); echo "   CONFLICT($rc) $f"
  fi
  cp "$W/out" "$f"
done
echo "merged: $n_clean clean, $n_conf with conflicts"
# step 1: the lane's new files, its current versions
for f in R/se-check.R tests/testthat/test-se-check.R \
         tests/testthat/test-allfit-start.R; do
  tr -d '\r' < "$LANE/$f" > "$f"; echo "   copied $f"
done
# steps 3 to 5 and 8 to 14
"$RS" dev/rel068-nanse-resolve.R || exit 1
grep -rln '^<<<<<<< \|^>>>>>>> ' $files R/se-check.R && {
  echo "CONFLICT MARKERS LEFT"; exit 1; }
# steps 2 and 7: NEWS, with nanse's bullets first
"$RS" dev/rel068-news-merge.R || exit 1
# the lane's untracked evidence (dev/), without overwriting
lst=$W/untracked.txt; : > "$lst"; clash=0
while IFS= read -r f; do
  case "$f" in dev/*) ;; *) continue ;; esac
  case "$f" in *stan-cache*/*) continue ;; esac
  if [ -e "$f" ]; then clash=$((clash + 1)); continue; fi
  printf '%s\n' "$f" >> "$lst"
done < <(git -C "$LANE" ls-files --others --exclude-standard)
( cd "$LANE" && tar -cf - -T "$lst" ) | tar -xpf -
echo "untracked dev files copied: $(grep -c '' "$lst"), already there: $clash"
echo "NANSE MERGE DONE"

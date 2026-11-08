#!/bin/bash
# Review of lane nlstart: the extension copies of the census test, against
# HEAD byte for byte, and their line endings.
cd /c/Users/adf44/source/r/frmtmb-wt-nlstart
export GIT_OPTIONAL_LOCKS=0
for e in coupling eam latent learn ode sample survival spatial; do
  f=extensions/frmtmb.$e/tests/testthat/test-conditions-census.R
  if [ ! -f "$f" ]; then echo "$e: no copy"; continue; fi
  h=dev/review-nlstart-out/census-head-$e.R
  git show "HEAD:$f" > "$h"
  if cmp -s "$f" "$h"; then s=same; else s=DIFF; fi
  cr=$(tr -cd '\r' < "$f" | wc -c)
  hcr=$(tr -cd '\r' < "$h" | wc -c)
  echo "$e: worktree vs HEAD bytes $s; CR in worktree $cr, in HEAD $hcr"
done
echo "git diff of extensions (content):"
git diff --stat -- extensions/
f=tests/testthat/test-conditions-census.R
echo "core copy CR: $(tr -cd '\r' < $f | wc -c)"

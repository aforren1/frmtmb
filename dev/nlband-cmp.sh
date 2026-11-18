#!/usr/bin/env bash
# Lane nlband: compare the ROW lines of two runs of dev/nlband-guard.R
# tag by tag: SAME where the fit's printed values (logLik and estimates to
# 17 digits, standard errors to 12, condition classes) are identical,
# else both lines. The HOLES count is the after arm's.
# Usage: bash dev/nlband-cmp.sh <before log> <after log>
b="$1"; a="$2"
grep "^ROW" "$a" | while IFS= read -r la; do
  tag="${la%% | *}"
  lb=$(grep -F "$tag |" "$b" | head -1)
  va="${la#* | }"; va="${va%% | HOLES*}"
  vb="${lb#* | }"; vb="${vb%% | HOLES*}"
  h="${la##*| HOLES }"
  if [ "$va" = "$vb" ]; then
    echo "SAME  $tag | holes $h"
  else
    echo "DIFF  $tag | holes $h"
    echo "  before: $vb" | cut -c1-200
    echo "  after:  $va" | cut -c1-200
  fi
done

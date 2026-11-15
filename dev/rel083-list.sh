#!/usr/bin/env bash
# The paths a lane tree holds beyond its HEAD, read without git status or
# git diff (both refresh the index's stat data): git ls-files lists
# modified, deleted and untracked-not-ignored paths by stat, and a
# modified path is kept only when its content (git hash-object, which
# applies the path's CRLF filter, no -w) differs from the index blob.
# Prints "<status> <path>" with status ?? (untracked), D (deleted) or M.
#
#   bash dev/rel083-list.sh <tree>
set -euo pipefail
export GIT_OPTIONAL_LOCKS=0
cd "$1"
git ls-files -z --others --exclude-standard | tr '\0' '\n' | sed '/^$/d' | sed 's/^/?? /'
del=$(git ls-files -z --deleted | tr '\0' '\n' | sed '/^$/d' | sort -u)
[ -z "$del" ] || printf '%s\n' "$del" | sed 's/^/D  /'
mod=$(git ls-files -z --modified | tr '\0' '\n' | sed '/^$/d' | sort -u)
if [ -n "$del" ]; then mod=$(comm -23 <(printf '%s\n' "$mod") <(printf '%s\n' "$del")); fi
[ -n "$mod" ] || exit 0
paste -d '\t' <(printf '%s\n' "$mod") <(printf '%s\n' "$mod" | git hash-object --stdin-paths) |
  awk -F'\t' 'NR == FNR { split($1, m, " "); blob[$2] = m[2]; next }
              !($1 in blob) || blob[$1] != $2 { print " M " $1 }' \
    <(git ls-files -s -z | tr '\0' '\n') -

#!/usr/bin/env bash
# Fail a job whose tests skipped because a Suggests package is broken
# or missing. testthat's skip_if_not_installed() turns both into a skip,
# and a skip passes R CMD check and covr, so a runner without a system
# library (fmesher without udunits2, 0.74.0) or an install.packages()
# that only warned hid real tests behind a green job.
#
#   skip-gate.sh <dir> [package allowed to be absent ...]
#
# <dir> is searched for testthat.Rout and testthat.Rout.fail (an
# R CMD check directory, or covr's install_path). A "{pkg} cannot be
# loaded" skip always fails: the package is installed and broken. A
# "{pkg} is not installed" skip fails unless pkg is named after <dir>,
# for a job that leaves a Suggests out on purpose.
set -u
dir=$1
shift
allowed=" $* "

# a read loop, not mapfile: core's R-CMD-check job runs this on macOS,
# whose own /bin/bash is 3.2 and has no mapfile
outs=()
while IFS= read -r f; do
  outs+=("$f")
done < <(find "$dir" -name 'testthat.Rout*' -type f 2>/dev/null)
if [ ${#outs[@]} -eq 0 ]; then
  echo "skip-gate: no testthat.Rout under $dir, so nothing was checked." >&2
  exit 1
fi
# An empty or cut-off output has no skip lines either, and read as
# clean: testthat's summary line must be there
summary='\[ FAIL [0-9]+ \| WARN [0-9]+ \| SKIP [0-9]+ \| PASS [0-9]+ \]'
if ! grep -qE "$summary" "${outs[@]}"; then
  echo "::error title=No testthat summary::no testthat summary line" \
    "([ FAIL n | WARN n | SKIP n | PASS n ]) in ${outs[*]}, so the" \
    "output is empty or cut off and nothing was checked."
  exit 1
fi

# testthat lists each skip reason once, as "• {pkg} <reason> (n): files";
# the bullet is "*" in a non-UTF-8 locale
lines=$(grep -hE '^(•|\*) \{[^}]+\} (cannot be loaded|is not installed)' \
  "${outs[@]}" | sort -u)

bad=0
while IFS= read -r l; do
  [ -z "$l" ] && continue
  pkg=$(sed -E 's/^[^{]*\{([^}]+)\}.*/\1/' <<< "$l")
  if [[ $l == *"is not installed"* && $allowed == *" $pkg "* ]]; then
    echo "skip-gate: allowed absent: $l"
    continue
  fi
  echo "::error title=Test skipped on a Suggests package::$l"
  bad=1
done <<< "$lines"

if [ $bad -eq 0 ]; then
  echo "skip-gate: no skip on a broken or missing Suggests package" \
    "(${#outs[@]} testthat output file(s))."
fi
exit $bad

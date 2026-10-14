#!/usr/bin/env bash
# Run every frmtmb.eam test file against the worktree source, one file
# per R process, two processes at a time. Logs go to
# dev/famuni-out/suite-<mode>/<file>.txt; the RESULT line of each is the
# count. Usage: bash dev/famuni-wiener-suite.sh [worktree|rellib] [files]
mode="${1:-worktree}"
shift
wt=/c/Users/adf44/source/r/frmtmb-wt-famuni
out="$wt/dev/famuni-out/suite-$mode"
mkdir -p "$out"
if [ "$#" -gt 0 ]; then
  files="$*"
else
  files=$(ls "$wt"/extensions/frmtmb.eam/tests/testthat/test-*.R)
fi
run_one() {
  f="$1"
  b=$(basename "$f" .R)
  "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" "$wt/dev/famuni-wiener-runtest.R" \
    "$mode" "$f" > "$out/$b.txt" 2>&1
}
export -f run_one
export mode out wt
printf '%s\n' $files | xargs -P 2 -I{} bash -c 'run_one {}'
grep -h "^RESULT" "$out"/*.txt
echo "files: $(printf '%s\n' $files | wc -l)  results: $(grep -l '^RESULT' "$out"/*.txt | wc -l)"

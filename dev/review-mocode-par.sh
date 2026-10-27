#!/usr/bin/env bash
# Runs every file in dev/review-mocode-log/rv-{core,ext}-files.txt
# against the reviewer's prototype build, one file per R process, five
# at a time. Logs: dev/review-mocode-log/rv-files/.
cd /c/Users/adf44/source/r/frmtmb-wt-mocode || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
out=dev/review-mocode-log/rv-files
mkdir -p "$out"
one() {
  f="$1"
  case "$f" in
    extensions/*) pkg=$(echo "$f" | cut -d/ -f2) ;;
    *) pkg=frmtmb ;;
  esac
  log="$out/$pkg-$(basename "$f" .R).txt"
  "$R" dev/review-mocode-run-test.R "$pkg" "$f" > "$log" 2>&1
}
export -f one
export R out
cat dev/review-mocode-log/rv-core-files.txt \
    dev/review-mocode-log/rv-ext-files.txt | xargs -P 5 -I{} bash -c 'one "$@"' _ {}
echo "ran $(ls "$out" | wc -l)"
grep -h '^RESULT' "$out"/*.txt | awk '{split($3,a,"=");split($4,b,"=");if(a[2]>0||b[2]>0)print}'

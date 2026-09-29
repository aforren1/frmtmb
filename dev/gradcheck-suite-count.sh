#!/bin/sh
# Totals from the per-file suite logs, which PowerShell writes as UTF-16.
# The count comes from RESULT lines, and a file with no result line is
# named rather than silently dropped.
set -e
cd "$(dirname "$0")/.."
arm="${1:-lane}"
dir="dev/gradcheck-log/suite-$arm"
flat="dev/gradcheck-log/suite-$arm-flat.txt"
: > "$flat"
missing=0
nfile=0
for f in "$dir"/*.txt; do
  nfile=$((nfile + 1))
  line=$(iconv -f UTF-16LE -t UTF-8 "$f" 2>/dev/null | sed 's/\r//' \
    | grep -E 'pass=|ABORTED' || true)
  if [ -z "$line" ]; then
    echo "NO RESULT $(basename "$f")" >> "$flat"
    missing=$((missing + 1))
  else
    printf '%s\n' "$line" >> "$flat"
  fi
done
echo "files: $nfile   with no result: $missing"
echo "--- files that are not clean:"
grep -v 'fail=0 error=0' "$flat" || echo "(none)"
echo "--- files with skips:"
grep -v 'skip=0' "$flat" | grep 'pass=' || echo "(none)"
echo "--- totals:"
awk '{for(i=1;i<=NF;i++){
        if($i ~ /^pass=/){split($i,a,"=");p+=a[2]}
        if($i ~ /^fail=/){split($i,a,"=");f+=a[2]}
        if($i ~ /^error=/){split($i,a,"=");e+=a[2]}
        if($i ~ /^skip=/){split($i,a,"=");s+=a[2]}}}
      END{print "PASS",p," FAIL",f," ERROR",e," SKIP",s}' "$flat"

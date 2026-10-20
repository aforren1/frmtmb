#!/usr/bin/env bash
# Reviewer of lane ndaterm: run test files one per process, P at a time.
#   bash dev/ndaterm-rev-run.sh <tag> <P> <rscript> <arm> <pkg:file> ...
# Logs dev/ndaterm-rev-out/<tag>/<pkg>__<file>.txt, summary RESULTS.txt.
set -u
TAG="$1"; P="$2"; RS="$3"; ARM="$4"; shift 4
WT=/c/Users/adf44/source/r/frmtmb-wt-ndaterm
OUT="$WT/dev/ndaterm-rev-out/$TAG"
mkdir -p "$OUT"
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export NOT_CRAN=true OPENBLAS_NUM_THREADS=4
run_one() {
  pf="$1"; pkg="${pf%%:*}"; f="${pf#*:}"
  log="$OUT/${pkg}__$(basename "$f" .R).txt"
  "$RS" "$WT/dev/ndaterm-runfile.R" "$pkg" "$WT/$f" "$ARM" > "$log" 2>&1
  grep -h "^RESULT\|^lib:" "$log" | tr '\n' ' '; echo
}
export -f run_one; export OUT WT RS ARM
printf '%s\n' "$@" | xargs -P "$P" -I{} bash -c 'run_one "$@"' _ {} \
  > "$OUT/RESULTS.txt"
echo done >> "$OUT/RESULTS.txt"

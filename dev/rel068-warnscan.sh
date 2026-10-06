#!/usr/bin/env bash
# dev/warnleak-scan.R over every test file of all eight packages on the
# release tree, twenty at a time: ungated, then the gated list with its
# gates set. Writes dev/release/warnscan-068/escaped.tsv (one line per
# escaped warning) and scanned.tsv (tier, file, lines written; a file
# the scan could not run writes the one line NOFILE, and nofile.txt
# names it).
#
#   bash dev/rel068-warnscan.sh
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
S=$ROOT/dev/rel068-log
OUT=$ROOT/dev/release/warnscan-068
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export NOT_CRAN=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
rm -rf "$OUT"; mkdir -p "$OUT/gated"
cd "$ROOT" || exit 1
run() {
  local tier="$1" list="$2"
  while read -r p f; do
    while [ "$(jobs -rp | wc -l)" -ge 12 ]; do sleep 1; done
    o="$OUT/$tier/$p--$(basename "$f" .R).tsv"
    ( "$RS" dev/warnleak-scan.R "$p" "$f" "$o" >/dev/null 2>&1
      [ -f "$o" ] || echo "NOFILE" > "$o" ) &
  done < "$list"
  wait
}
# every file with every gate set: the gated run is a superset of the
# ungated one (0.68.0)
FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true FRMTMB_FUZZ=true \
  run gated "$S/suite.jobs"
: > "$OUT/scanned.tsv"; : > "$OUT/escaped.tsv"; : > "$OUT/nofile.txt"
for t in gated; do
  for f in "$OUT/$t"/*.tsv; do
    printf '%s\t%s\t%s\n' "$t" "$(basename "$f" .tsv)" \
      "$(grep -c '' "$f")" >> "$OUT/scanned.tsv"
    grep -q "^NOFILE$" "$f" &&
      echo "$t $(basename "$f" .tsv)" >> "$OUT/nofile.txt"
    grep -v '^NOFILE$' "$f" | grep -v '^$' >> "$OUT/escaped.tsv"
  done
done
rm -rf "$OUT/ungated" "$OUT/gated"
echo "scanned $(grep -c '' "$OUT/scanned.tsv") files," \
  "escaped $(grep -c '' "$OUT/escaped.tsv")"
echo "WARNSCAN DONE"

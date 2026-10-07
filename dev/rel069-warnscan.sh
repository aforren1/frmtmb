#!/usr/bin/env bash
# dev/warnleak-scan.R over every test file of all eight packages on the
# release tree, every gate set (brms, drmTMB, fuzz, scale), N at a time,
# with dev/rel069-profile.R as R_PROFILE_USER: one run gives the list of
# escaped warnings and the warning-precedence log (one report per
# parameter). Writes dev/release/warnscan-069/escaped.tsv (one line per
# escaped warning), scanned.tsv (file, lines written; a file the scan
# could not run writes the one line NOFILE, and nofile.txt names it) and
# dev/rel069-log/prof/<pid>.txt, summarized by dev/rel069-prof-sum.R.
#
#   bash dev/rel069-warnscan.sh <job list> [N]
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
list="$1"; N="${2:-16}"
OUT=$ROOT/dev/release/warnscan-069
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
RT=/c/rtools45
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true FRMTMB_DRMTMB_FIT_TESTS=true
export FRMTMB_FUZZ=true FRMTMB_SCALE_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export PATH="$RT/usr/bin:$RT/x86_64-w64-mingw32.static.posix/bin:$PATH"
export R_PROFILE_USER="C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel069-profile.R"
export REL069_PROF_LOG="C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel069-log/prof"
: "${TMP:?TMP is unset}"
rm -rf "$OUT" "$REL069_PROF_LOG"; mkdir -p "$OUT/files" "$REL069_PROF_LOG"
cd "$ROOT" || exit 1
while read -r p f; do
  [ -n "$p" ] || continue
  while [ "$(jobs -rp | wc -l)" -ge "$N" ]; do sleep 1; done
  o="$OUT/files/$p--$(basename "$f" .R).tsv"
  ( REL069_PROF_FILE="$p/$(basename "$f")" \
      "$RS" dev/warnleak-scan.R "$p" "$f" "$o" >/dev/null 2>&1
    [ -f "$o" ] || echo "NOFILE" > "$o" ) &
done < "$list"
wait
: > "$OUT/scanned.tsv"; : > "$OUT/escaped.tsv"; : > "$OUT/nofile.txt"
for f in "$OUT/files"/*.tsv; do
  printf '%s\t%s\n' "$(basename "$f" .tsv)" "$(grep -c '' "$f")" \
    >> "$OUT/scanned.tsv"
  grep -q "^NOFILE$" "$f" && basename "$f" .tsv >> "$OUT/nofile.txt"
  grep -v '^NOFILE$' "$f" | grep -v '^$' >> "$OUT/escaped.tsv"
done
echo "scanned $(grep -c '' "$OUT/scanned.tsv") of $(grep -c '' "$list")" \
  "files, escaped $(grep -c '' "$OUT/escaped.tsv")," \
  "could not run $(grep -c '' "$OUT/nofile.txt")"
echo "WARNSCAN DONE"

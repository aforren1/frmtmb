#!/usr/bin/env bash
# Review: a spot-check of the floors, by dev/rel076-floor.sh's method
# (load a COPY of an installed build; nothing is installed), but with
# every write kept to the review's own places: the copies go to the
# review library C:/Users/adf44/source/r/rev-rel076-lib (deleted at the
# end), the logs to dev/review-rel076-out/floor/. dev/rel076-floor.sh
# itself is not run: it writes dev/rel076-log/floor-*.jobs and, through
# dev/ciharden-scan.sh, dev/ciharden-log/r076-floor-*, the
# consolidation's own evidence.
#   bash dev/review-rel076-floor.sh
REL=/c/Users/adf44/source/r/frmtmb-wt-rel076
SCR=C:/Users/adf44/source/r/rev-rel076-lib
R13=C:/Users/adf44/source/r/rellib-r13
R14=C:/Users/adf44/source/r/rellib-r14
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
OUT=$REL/dev/review-rel076-out/floor
mkdir -p "$OUT" "$SCR"
export NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$REL/dev/review-rel076-out/stan-cache
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f)"
}
echo "== exports at 70476dc3 (any line of NAMESPACE, multi-line blocks included)"
for s in autocor_flat_priors autocor_bound_chart autocor_innov_sd autocor_row_loglik; do
  printf '  %s: 70476dc3 %s, tree %s\n' "$s" \
    "$(git -C "$REL" show 70476dc3:NAMESPACE | grep -cw "$s")" \
    "$(grep -cw "$s" "$REL/NAMESPACE")"
done
for p in frmtmb.spatial frmtmb.sample; do
  echo "== $p $(grep '^Version:' "$REL/extensions/$p/DESCRIPTION" | cut -d' ' -f2), $(grep -o 'frmtmb (>= [0-9.]*)' "$REL/extensions/$p/DESCRIPTION" | head -1)"
  rm -rf "${SCR:?}/old-$p"; mkdir -p "$SCR/old-$p"; cp -r "$R14/$p" "$SCR/old-$p/"
  lower_floor "$SCR/old-$p/$p"
  echo "-- old: rellib-r14's build over rellib-r13, floor lowered"
  out=$(load_line "$SCR/old-$p" "$R13" "$p"); echo "$out"
  if echo "$out" | grep -q LOADED && [ "$p" = frmtmb.sample ]; then
    for f in test-default-priors-autocor test-loglik-autocor test-mo-coding test-loo; do
      ROUND_BASE=$R13 "$RS" "$REL/dev/ciharden-run1.R" "$SCR/old-$p" "$p" \
        "$REL/extensions/$p/tests/testthat/$f.R" > "$OUT/old-$p-$f.log" 2>&1
      echo "  $(grep -a '^lib:' "$OUT/old-$p-$f.log") | $(grep -a '^RESULT' "$OUT/old-$p-$f.log")"
      grep -a "could not find function" "$OUT/old-$p-$f.log" | sed 's/.*could not find function/  could not find function/' | sort | uniq -c | head -5
    done
  fi
  echo "-- new: rellib-r14 as built"
  load_line "$R14" "$R14" "$p"
  rm -rf "${SCR:?}/conv-$p"; mkdir -p "$SCR/conv-$p"; cp -r "$R13/$p" "$SCR/conv-$p/"
  echo "-- converse: $p $(grep '^Version:' "$SCR/conv-$p/$p/DESCRIPTION" | cut -d' ' -f2) from rellib-r13 over rellib-r14"
  load_line "$SCR/conv-$p" "$R14" "$p"
done
rm -rf C:/Users/adf44/source/r/rev-rel076-lib
echo "FLOOR SPOT-CHECK DONE"

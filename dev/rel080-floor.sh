#!/usr/bin/env bash
# The floors of the 0.80.0 release, measured, by dev/rel079-floor.sh's
# method: nothing is installed; each arm loads a COPY of an installed build
# in front of the core it is measured on.
#   old: rellib-r19's frmtmb.sample and frmtmb.latent, each floor lowered
#        to 0.61.0, over rellib-r18 (0.79.1); every test file of the suite
#        from the release tree, one per process, every gate on.
#   new: rellib-r19 as built. converse: rellib-r18's build over rellib-r19.
# The six other extensions: every test file from the release tree on
# rellib-r18 alone, every gate on.
#
#   bash dev/rel080-floor.sh [N] [part]   part: frmtmb.sample, frmtmb.latent, others
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/2}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/2}"
export GIT_OPTIONAL_LOCKS=0
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel080
SCR=C:/Users/adf44/source/r/rel080-floor-lib
RBASE=C:/Users/adf44/source/r/rellib-r18
RNEW=C:/Users/adf44/source/r/rellib-r19
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
N="${1:-8}"
part="${2:-all}"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "== 3cd2a0c8 (0.79.1): the seams the new floors rest on"
printf "  logLik_penalized in NAMESPACE: 3cd2a0c8 %s, this tree %s\n" \
  "$(git show 3cd2a0c8:NAMESPACE | grep -c logLik_penalized)" \
  "$(grep -c logLik_penalized NAMESPACE)"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f); cat('  floor in Meta/package.rds:', paste(vapply(c(m\$Depends, m\$Imports), function(d) paste(d\$name, format(d\$version)), ''), collapse = ', '), '\n')"
}
for p in frmtmb.sample frmtmb.latent; do
  [ "$part" = all ] || [ "$part" = "$p" ] || continue
  v=$(grep '^Version:' "extensions/$p/DESCRIPTION" | cut -d' ' -f2 | tr -d '\r')
  fl=$(grep -o 'frmtmb (>= [0-9.]*)' "extensions/$p/DESCRIPTION" | head -1)
  echo "== $p $v, floor $fl"
  FL=$SCR/old-$p; rm -rf "${SCR:?}/old-$p"; mkdir -p "$FL"
  cp -r "$RNEW/$p" "$FL/"
  echo "-- old: rellib-r19's build over rellib-r18 (floor lowered)"
  lower_floor "$FL/$p"
  out=$(load_line "$FL" "$RBASE" "$p"); echo "$out"
  if echo "$out" | grep -q "LOADED"; then
    : > "dev/rel080-log/floor-$p.jobs"
    for f in extensions/$p/tests/testthat/test-*.R; do
      echo "$p $f" >> "dev/rel080-log/floor-$p.jobs"
    done
    ROUND_BASE=$RBASE /usr/bin/bash dev/ciharden-scan.sh \
      "r080-floor-$p" ref "$FL" "$N" "dev/rel080-log/floor-$p.jobs"
  fi
  echo "-- new: rellib-r19 as built"
  load_line "$RNEW" "$RNEW" "$p"
  FC=$SCR/conv-$p; rm -rf "${SCR:?}/conv-$p"; mkdir -p "$FC"
  cp -r "$RBASE/$p" "$FC/"
  echo "-- converse: $p $(grep '^Version:' "$FC/$p/DESCRIPTION" | cut -d' ' -f2) from rellib-r18 over rellib-r19"
  load_line "$FC" "$RNEW" "$p"
done
if [ "$part" = all ] || [ "$part" = others ]; then
  echo "== the six other extensions: the release tree's tests on rellib-r18"
  : > dev/rel080-log/floor-others.jobs
  for e in coupling eam learn ode spatial survival; do
    for f in extensions/frmtmb.$e/tests/testthat/test-*.R; do
      echo "frmtmb.$e $f" >> dev/rel080-log/floor-others.jobs
    done
  done
  ROUND_BASE=$RBASE /usr/bin/bash dev/ciharden-scan.sh r080-floor-others ref \
    base "$N" dev/rel080-log/floor-others.jobs
fi
rm -rf "$SCR"
echo "FLOOR DONE"

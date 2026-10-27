#!/usr/bin/env bash
# The floors of the 0.76.0 release, measured. Adapted from
# dev/rel075-floor.sh, with one change of method: nothing is installed.
# The consolidation installs only into rellib-r14, so each arm loads a
# COPY of an installed build (none of these packages has compiled code)
# in front of the core it is measured on.
#
# Arm "old": rellib-r14's build of each package whose floor is raised to
# frmtmb (>= 0.76.0), copied into a scratch library with its declared
# floor lowered to 0.61.0 (DESCRIPTION and Meta/package.rds, so that
# library() does not refuse on the version), loaded in front of
# rellib-r13 (frmtmb 0.75.0 and the 0.75.0 extensions). frmtmb.spatial
# registers sar() with frmtmb_register_autocor(jacobian =), which 0.75.0
# does not take, so it must fail to load there. frmtmb.sample: if it
# loads, every test file of its suite runs there (the release tree's
# files), one per process, every gate on; a file that fails only on the
# old core names the floor.
# Arm "new": rellib-r14 as built, loaded.
# Arm "converse": rellib-r13's 0.75.0 build of each of these packages
# (frmtmb.sample 0.22.0, frmtmb.spatial 0.2.0) loaded in front of
# rellib-r14's frmtmb 0.76.0, to say whether an old extension still
# loads on the new core.
#
#   bash dev/rel076-floor.sh [N]
#
# The scratch libraries live outside the tree, in
# C:/Users/adf44/source/r/rel076-floor-lib, and are removed at the end.
# PATH is left as Git's bash gives it (dev/rel075-floor.sh says why).
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/12}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/12}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel076
SCR=C:/Users/adf44/source/r/rel076-floor-lib
R13=C:/Users/adf44/source/r/rellib-r13
R14=C:/Users/adf44/source/r/rellib-r14
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "PATH starts: $(echo "$PATH" | cut -d: -f1-3); bash: $BASH"
echo "== 70476dc3 (0.75.0): the seams the new floors rest on"
for s in autocor_flat_priors autocor_bound_chart autocor_innov_sd autocor_row_loglik; do
  printf '  NAMESPACE export(%s): 70476dc3 %s, this tree %s\n' "$s" \
    "$(git show 70476dc3:NAMESPACE | grep -c "^export($s)")" \
    "$(grep -c "^export($s)" NAMESPACE)"
done
printf '  jacobian = in frmtmb_register_autocor (R/autocor-seam.R): 70476dc3 %s, this tree %s\n' \
  "$(git show 70476dc3:R/autocor-seam.R | grep -c 'jacobian = NULL')" \
  "$(grep -c 'jacobian = NULL' R/autocor-seam.R)"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
lower_floor() {   # <installed package dir>
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$1/DESCRIPTION"
  "$RS" -e "f <- '$1/Meta/package.rds'; m <- readRDS(f); for (k in c('Depends', 'Imports')) m[[k]] <- lapply(m[[k]], function(d) { if (identical(d\$name, 'frmtmb') && !is.null(d\$version)) d\$version <- package_version('0.61.0'); d }); saveRDS(m, f); cat('  floor in Meta/package.rds:', paste(vapply(m\$Depends, function(d) paste(d\$name, format(d\$version)), ''), collapse = ', '), '\n')"
}
for p in frmtmb.spatial frmtmb.sample; do
  v=$(grep '^Version:' "extensions/$p/DESCRIPTION" | cut -d' ' -f2)
  fl=$(grep -o 'frmtmb (>= [0-9.]*)' "extensions/$p/DESCRIPTION" | head -1)
  echo "== $p $v, floor $fl"
  FL=$SCR/old-$p; rm -rf "${SCR:?}/old-$p"; mkdir -p "$FL"
  cp -r "$R14/$p" "$FL/"
  echo "-- old: rellib-r14's build over rellib-r13 (floor lowered)"
  lower_floor "$FL/$p"
  out=$(load_line "$FL" "$R13" "$p"); echo "$out"
  if echo "$out" | grep -q "LOADED" && [ "$p" = frmtmb.sample ]; then
    : > "dev/rel076-log/floor-$p.jobs"
    for f in extensions/$p/tests/testthat/test-*.R; do
      echo "$p $f" >> "dev/rel076-log/floor-$p.jobs"
    done
    ROUND_BASE=$R13 /usr/bin/bash dev/ciharden-scan.sh \
      "r076-floor-$p" ref "$FL" "${1:-8}" "dev/rel076-log/floor-$p.jobs"
  fi
  echo "-- new: rellib-r14 as built"
  load_line "$R14" "$R14" "$p"
  FC=$SCR/conv-$p; rm -rf "${SCR:?}/conv-$p"; mkdir -p "$FC"
  cp -r "$R13/$p" "$FC/"
  echo "-- converse: $p $(grep '^Version:' "$FC/$p/DESCRIPTION" | cut -d' ' -f2) from rellib-r13 over rellib-r14"
  load_line "$FC" "$R14" "$p"
done
rm -rf C:/Users/adf44/source/r/rel076-floor-lib
echo "FLOOR DONE"

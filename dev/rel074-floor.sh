#!/usr/bin/env bash
# The floors of the 0.74.0 release, measured. Adapted from
# dev/rel073-floor.sh, with an arm per package whose floor is claimed.
#
# Arm "old": each package is installed ALONE into a scratch library over
# rellib-r11 (frmtmb 0.73.0 and the 0.73.0 extensions), from a copy whose
# DESCRIPTION floor is lowered so that R does not refuse on the version,
# then loaded. frmtmb.survival, frmtmb.spatial and frmtmb.latent claim
# frmtmb (>= 0.74.0), so each must fail there. frmtmb.sample and
# frmtmb.learn: if one installs and loads, every test file of its suite
# runs there, one per process, every gate on; a file that fails only on
# the old core names the floor.
# Arm "new": the same packages, unmodified, installed alone into a
# scratch library over rellib-r12 (frmtmb 0.74.0), then loaded.
# frmtmb 0.73.1 (19b162b9) was never installed as a library; it differs
# from 0.73.0 only in R/predict.R, and its NAMESPACE and frmtmb_family()
# are read from git below to show it lacks the same seams.
#
#   bash dev/rel074-floor.sh [N]
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/3}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/3}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-release
FLD=C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel074-log/floor-lib
R11=C:/Users/adf44/source/r/rellib-r11
R12=C:/Users/adf44/source/r/rellib-r12
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
cd "$ROOT" || exit 1
echo "== 19b162b9 (0.73.1): the seams the new floors rest on"
for s in frmtmb_register_covstruct us_sigma us_theta frmtmb_local_registries; do
  printf '  NAMESPACE export(%s): %s\n' "$s" \
    "$(git show 19b162b9:NAMESPACE | grep -c "^export($s)")"
done
printf '  frmtmb_family( simplex_pars = / aterm_data = in R/families.R: %s / %s\n' \
  "$(git show 19b162b9:R/families.R | grep -c 'simplex_pars = ')" \
  "$(git show 19b162b9:R/families.R | grep -c 'aterm_data = NULL')"
printf '  this tree: %s / %s\n' "$(grep -c 'simplex_pars = ' R/families.R)" \
  "$(grep -c 'aterm_data = NULL' R/families.R)"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
inst_line() {   # <lib> <base> <src>
  R_LIBS="$1;$2;$USR" "$R" CMD INSTALL --library="$1" --no-multiarch "$3" 2>&1 |
    grep -E 'DONE|ERROR|Error|error|object .* not exported|is not an exported|installation of package|unused argument' |
    sed 's/^/  /' | head -8
}
for p in frmtmb.survival frmtmb.spatial frmtmb.latent frmtmb.sample frmtmb.learn; do
  v=$(grep '^Version:' "extensions/$p/DESCRIPTION" | cut -d' ' -f2)
  fl=$(grep -o 'frmtmb (>= [0-9.]*)' "extensions/$p/DESCRIPTION" | head -1)
  echo "== $p $v, floor $fl"
  # arm old
  FL=$FLD/old-$p; rm -rf "$FL"; mkdir -p "$FL"
  SRC=$FLD-src/$p; rm -rf "$SRC"; mkdir -p "$FLD-src"
  cp -r "extensions/$p" "$SRC"
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$SRC/DESCRIPTION"
  echo "-- old: install over rellib-r11 (floor lowered: $(grep -o 'frmtmb (>= [0-9.]*)' "$SRC/DESCRIPTION" | head -1))"
  inst_line "$FL" "$R11" "$SRC"
  if [ -d "$FL/$p" ]; then
    load_line "$FL" "$R11" "$p"
    case $p in frmtmb.sample|frmtmb.learn)
      : > "dev/rel074-log/floor-$p.jobs"
      for f in extensions/$p/tests/testthat/test-*.R; do
        echo "$p $f" >> "dev/rel074-log/floor-$p.jobs"
      done
      ROUND_BASE=$R11 /usr/bin/bash dev/ciharden-scan.sh \
        "r074-floor-$p" ref "$FL" "${1:-8}" "dev/rel074-log/floor-$p.jobs";;
    esac
  else
    echo "  NOT INSTALLED: $p $v does not install over frmtmb 0.73.0"
  fi
  # arm new
  FN=$FLD/new-$p; rm -rf "$FN"; mkdir -p "$FN"
  echo "-- new: install over rellib-r12 (floor as written)"
  inst_line "$FN" "$R12" "extensions/$p"
  load_line "$FN" "$R12" "$p"
done
echo "FLOOR DONE"

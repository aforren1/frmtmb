#!/usr/bin/env bash
# The floors of the 0.75.0 release, measured. Adapted from
# dev/rel074-floor.sh.
#
# Arm "old": each package whose floor is raised to frmtmb (>= 0.75.0) is
# installed ALONE into a scratch library over rellib-r12 (frmtmb 0.74.0
# and the 0.74.0 extensions), from a copy whose DESCRIPTION floor is
# lowered so that R does not refuse on the version, then loaded.
# frmtmb.survival (frmtmb_register_aterm(option =) in .onLoad()) and
# frmtmb.spatial (it imports frmtmb_register_autocor()) must fail there.
# frmtmb.sample and frmtmb.latent: if one installs and loads, every test
# file of its suite runs there, one per process, every gate on; a file
# that fails only on the old core names the floor.
# Arm "new": the same packages, unmodified, installed alone into a
# scratch library over rellib-r13 (frmtmb 0.75.0), then loaded.
# Arm "converse": each 0.74.0 build of these four (from rellib-r12) is
# loaded in front of rellib-r13's frmtmb 0.75.0, to say whether an old
# extension still loads on the new core.
# frmtmb 0.74.0 is 18b536fc's and 4fd134db's code (the CI commit after
# it changed no package file), so rellib-r12 is the old core.
#
#   bash dev/rel075-floor.sh [N]
#
# The scratch libraries and the source copies live outside the tree,
# in C:/Users/adf44/source/r/rel075-floor-lib, and are removed at the
# end; nothing is installed into the tree or into a release library.
# PATH is left as Git's bash gives it: Rtools first on PATH made an R
# child segfault at exit in this round (dev/round-20261024.md, "The
# ported brms suite"), and the floor packages have no compiled code.
# The suite arm calls dev/ciharden-scan.sh, the harness every per-file
# run of the round used, which sets its own PATH inside its process.
set -u
export TMP="${TMP:-C:/Users/adf44/AppData/Local/Temp/12}" TEMP="${TEMP:-C:/Users/adf44/AppData/Local/Temp/12}"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel075
SCR=C:/Users/adf44/source/r/rel075-floor-lib
FLD=$SCR/lib
R12=C:/Users/adf44/source/r/rellib-r12
R13=C:/Users/adf44/source/r/rellib-r13
USR=C:/Users/adf44/AppData/Local/R/win-library/4.6
R="/c/Program Files/R/R-4.6.1/bin/R.exe"
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
cd "$ROOT" || exit 1
mkdir -p "$SCR"
echo "PATH starts: $(echo "$PATH" | cut -d: -f1-3); bash: $BASH"
echo "== 18b536fc (0.74.0): the seams the new floors rest on"
for s in frmtmb_register_autocor frmtmb_autocor_pars frmtmb_covstruct_info; do
  printf '  NAMESPACE export(%s): 18b536fc %s, this tree %s\n' "$s" \
    "$(git show 18b536fc:NAMESPACE | grep -c "^export($s)")" \
    "$(grep -c "^export($s)" NAMESPACE)"
done
printf '  frmtmb_register_aterm( option = in R/parse.R: 18b536fc %s, this tree %s\n' \
  "$(git show 18b536fc:R/parse.R | grep -c 'option = NULL')" \
  "$(grep -c 'option = NULL' R/parse.R)"
load_line() {   # <lib> <base> <pkg>
  "$RS" -e ".libPaths(c('$1', '$2', '$USR')); cat('  frmtmb', format(packageVersion('frmtmb')), dirname(find.package('frmtmb')), '\n'); r <- tryCatch({ suppressPackageStartupMessages(library('$3', character.only = TRUE)); paste('  LOADED', '$3', format(packageVersion('$3')), dirname(find.package('$3'))) }, error = function(e) paste('  LOAD ERROR:', conditionMessage(e))); cat(r, '\n')" 2>&1
}
inst_line() {   # <lib> <base> <src>
  R_LIBS="$1;$2;$USR" "$R" CMD INSTALL --library="$1" --no-multiarch "$3" 2>&1 |
    grep -E 'DONE|ERROR|Error|error|object .* not exported|is not an exported|installation of package|unused argument' |
    sed 's/^/  /' | head -8
}
for p in frmtmb.survival frmtmb.spatial frmtmb.sample frmtmb.latent; do
  v=$(grep '^Version:' "extensions/$p/DESCRIPTION" | cut -d' ' -f2)
  fl=$(grep -o 'frmtmb (>= [0-9.]*)' "extensions/$p/DESCRIPTION" | head -1)
  echo "== $p $v, floor $fl"
  # arm old
  FL=$FLD/old-$p; rm -rf "${FLD:?}/old-$p"; mkdir -p "$FL"
  SRC=$SCR/src/$p; rm -rf "${SCR:?}/src/$p"; mkdir -p "$SCR/src"
  cp -r "extensions/$p" "$SRC"
  sed -i -E 's/frmtmb \(>= [0-9.]+\)/frmtmb (>= 0.61.0)/' "$SRC/DESCRIPTION"
  echo "-- old: install over rellib-r12 (floor lowered: $(grep -o 'frmtmb (>= [0-9.]*)' "$SRC/DESCRIPTION" | head -1))"
  inst_line "$FL" "$R12" "$SRC"
  if [ -d "$FL/$p" ]; then
    load_line "$FL" "$R12" "$p"
    case $p in frmtmb.sample|frmtmb.latent)
      : > "dev/rel075-log/floor-$p.jobs"
      for f in extensions/$p/tests/testthat/test-*.R; do
        echo "$p $f" >> "dev/rel075-log/floor-$p.jobs"
      done
      ROUND_BASE=$R12 /usr/bin/bash dev/ciharden-scan.sh \
        "r075-floor-$p" ref "$FL" "${1:-8}" "dev/rel075-log/floor-$p.jobs";;
    esac
  else
    echo "  NOT INSTALLED: $p $v does not install over frmtmb 0.74.0"
  fi
  # arm new
  FN=$FLD/new-$p; rm -rf "${FLD:?}/new-$p"; mkdir -p "$FN"
  echo "-- new: install over rellib-r13 (floor as written)"
  inst_line "$FN" "$R13" "extensions/$p"
  load_line "$FN" "$R13" "$p"
  # arm converse: the 0.74.0 build in front of frmtmb 0.75.0
  FC=$FLD/conv-$p; rm -rf "${FLD:?}/conv-$p"; mkdir -p "$FC"
  cp -r "$R12/$p" "$FC/"
  echo "-- converse: $p $(grep '^Version:' "$FC/$p/DESCRIPTION" | cut -d' ' -f2) from rellib-r12 over rellib-r13"
  load_line "$FC" "$R13" "$p"
done
rm -rf C:/Users/adf44/source/r/rel075-floor-lib
echo "FLOOR DONE"

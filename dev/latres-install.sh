#!/usr/bin/env bash
# Roxygenise and install core (and any extension named) into the lane's
# private library. Usage: bash dev/latres-install.sh [--no-rox] [pkg ...]
source "$(dirname "$0")/latres-env.sh"
cd "$LATRES_ROOT" || exit 1
rox=1
if [ "${1:-}" = "--no-rox" ]; then rox=0; shift; fi
libs="c('$LATRES_LIB','$BASE_LIB','$USER_LIB')"
if [ $rox = 1 ]; then
  "$RSCRIPT" -e ".libPaths($libs); roxygen2::roxygenise('$LATRES_ROOT')" \
    > dev/latres-out/rox.log 2>&1 || { tail -20 dev/latres-out/rox.log; exit 1; }
fi
R_LIBS="$LATRES_LIB;$BASE_LIB;$USER_LIB" "$RCMD" CMD INSTALL \
  --library="$LATRES_LIB" --no-multiarch --no-docs "$LATRES_ROOT" \
  > dev/latres-out/install-core.log 2>&1 || {
    tail -20 dev/latres-out/install-core.log; exit 1; }
for p in "$@"; do
  R_LIBS="$LATRES_LIB;$BASE_LIB;$USER_LIB" "$RCMD" CMD INSTALL \
    --library="$LATRES_LIB" --no-multiarch --no-docs \
    "$LATRES_ROOT/extensions/$p" > "dev/latres-out/install-$p.log" 2>&1 || {
      tail -20 "dev/latres-out/install-$p.log"; exit 1; }
done
echo installed

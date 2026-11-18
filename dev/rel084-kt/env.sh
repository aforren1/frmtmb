# Environment for the rel084 key-table reruns: rellib-r23 build, release tree.
# Usage: source dev/rel084-kt/env.sh; then krun <log name> <script> [args]
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export TMP=C:/Users/adf44/AppData/Local/Temp/2
export TEMP=C:/Users/adf44/AppData/Local/Temp/2
export NOT_CRAN=true
export FRMTMB_STAN_CACHE=C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/stan-cache
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export R_LIBS=C:/Users/adf44/source/r/rellib-r23
export R_LIBS_USER=C:/Users/adf44/AppData/Local/R/win-library/4.6
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
ROOT=/c/Users/adf44/source/r/frmtmb-wt-rel084
KLOG=$ROOT/dev/rel084-log/kt
krun() {
  local log=$1; shift
  cd "$ROOT" || return 1
  echo "START $(date '+%F %T') $*" > "$KLOG/$log"
  "$RS" "$@" >> "$KLOG/$log" 2>&1
  local rc=$?
  echo "EXIT $rc $(date '+%F %T')" >> "$KLOG/$log"
  return $rc
}

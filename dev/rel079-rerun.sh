#!/usr/bin/env bash
# Each lane's key tables on the merged build. Every script runs from a copy
# under dev/rel079-log/rr/<lane>/, with the lane's library and the reviewer's
# library read as rellib-r17 and the lane's sourced helpers read from the same
# folder, so nothing is written in a lane's tree. A log's `lib:` line, written by
# the runner, names the build that frmtmb loads from.
#
#   bash dev/rel079-rerun.sh <lane> <script> [args...]
#
# One process. xargs -P runs several; the job list is dev/rel079-rerun.jobs.
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel079/dev/rel079-env.sh
ROOT=C:/Users/adf44/source/r/frmtmb-wt-rel079
cd "$ROOT" || exit 1
RS="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export R_MAKEVARS_USER=C:/Users/adf44/Documents/.R/Makevars.win
export FRMTMB_STAN_CACHE=$ROOT/dev/stan-cache
export NOT_CRAN=true
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
l="$1"; s="$2"; shift 2
tag=$(printf '%s' "$s $*" | sed 's#[ /:]#_#g' | cut -c1-120)
J=dev/rel079-log/rr/$l
D=$J/job-$(printf '%s' "$tag" | md5sum | cut -c1-12)
mkdir -p "$D/out"
LANE=C:/Users/adf44/source/r/frmtmb-wt-$l
ed() {
  sed -e "s#C:/Users/adf44/source/r/wt-[a-z0-9]*-lib[A-Za-z0-9/_-]*#C:/Users/adf44/source/r/rellib-r17#g" \
      -e "s#C:/Users/adf44/source/r/rev-$l-lib[A-Za-z0-9/_-]*#C:/Users/adf44/source/r/rellib-r17#g" \
      -e "s#source(\"dev/\(review-[a-z0-9-]*\.R\)\"#source(\"$D/\1\"#g" \
      -e "s#source(\"C:/Users/adf44/source/r/frmtmb-wt-$l/dev/\(review-[a-z0-9-]*\.R\)\"#source(\"$ROOT/$D/\1\"#g" \
      -e "s#paste0(root, \"\(review-$l-common\.R\)\")#\"$ROOT/$D/\1\"#g" \
      -e "s#frmtmb-wt-$l\"#frmtmb-wt-rel079\"#g" \
      -e "s#frmtmb-wt-$l/dev/stan-cache#frmtmb-wt-rel079/dev/stan-cache#g" \
      -e "s#dev/review-$l-out/stan-cache#dev/stan-cache#g" \
      -e "s#\"dev/review-$l-out/[A-Za-z0-9/_-]*\"#\"$D/out\"#g" "$1"
}
[ "$l" = rf1 ] && export REV_LIB=C:/Users/adf44/source/r/rellib-r17
export REV_R1_LIB=C:/Users/adf44/source/r/rellib-r17
for f in dev/review-$l-common.R dev/review-$l-lib.R dev/review-$l-designs*.R; do
  [ -f "$f" ] && ed "$f" > "$D/$(basename "$f")"
done
{
  # the process that ran the script says which build it read, at its end
  cat <<'RPRE'
reg.finalizer(globalenv(), function(e) {
  if ("frmtmb" %in% loadedNamespaces()) {
    p <- getNamespaceInfo("frmtmb", "path")
    cat("lib-used:", p, packageDescription("frmtmb", lib.loc = dirname(p))$Built,
        format(file.mtime(file.path(p, "DESCRIPTION"))), "\n")
  }
}, onexit = TRUE)
RPRE
  ed "dev/$s"
} > "$D/$s"
{
  echo "lib: reading $D/$s"
  arm=after; [ "$l" = latres ] && arm=lane
  REV_ARM="$arm" "$RS" "$D/$s" "$@"
  echo "EXIT $?"
} > "$J/$tag.log" 2>&1

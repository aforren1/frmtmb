#!/usr/bin/env bash
# Lane nanse's attack set (as lane setier ran it, dev/setier-p1-attacks.sh)
# on the 0.69.0 release library and on lane setier's final build, with
# the reference BLAS. Each script's "lane" arm is pointed at the library
# under test by a copy in the session scratchpad, so the scripts in dev/
# stay as the lane left them.
#
#   bash dev/rel069-attacks.sh
set -u
REL=/c/Users/adf44/source/r/frmtmb-wt-release
SP=/c/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/attacks
cd "$REL"
L=dev/rel069-log/attacks
mkdir -p "$L" "$SP"
REF="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
for s in dev/setier-rev-a-rev-cases.R dev/setier-rev-a-rev-spread2.R \
         dev/setier-rev-a-rev-grby.R dev/setier-rev-a-rev-pred.R \
         dev/setier-rev-a-rev2-b1.R dev/setier-rev-a-rev2-b1-rho.R \
         dev/setier-rev-a-rev3-falseloss-sweep.R \
         dev/setier-rev-a-rev2-ridge-re.R; do
  b=$(basename "$s" .R)
  sed 's#C:/Users/adf44/source/r/wt-setier-lib#C:/Users/adf44/source/r/rellib-r7#g' \
    "$s" > "$SP/$b-rel.R"
  "$REF" "$SP/$b-rel.R" lane > "$L/$b-rel-ref.txt" 2>&1 &
  "$REF" "$s" lane > "$L/$b-setier-ref.txt" 2>&1 &
done
wait
echo ATTACKS DONE

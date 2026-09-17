#!/usr/bin/env bash
# Reviewer, lane wt-conditions: run the lane's parse-tree proof against
# the BASE sources (the main checkout, at e031e8c) rather than the
# lane's own snapshot, then mutate one message and one argument in a
# copy of the lane tree and show the proof fails.
#   bash dev/conditions-rev-proof.sh
set -u
wt=/c/Users/adf44/source/r/frmtmb-wt-conditions
main=/c/Users/adf44/source/r/frmtmb
out=$wt/dev/conditions-rev-proof
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
rm -rf "$out"; mkdir -p "$out"
copy_tree() { # src dst
  mkdir -p "$2"
  cp -r "$1/R" "$1/inst" "$2/"
  for e in "$1"/extensions/*/; do
    n=$(basename "$e"); mkdir -p "$2/extensions/$n"
    cp -r "$e/R" "$2/extensions/$n/"
    [ -d "$e/inst" ] && cp -r "$e/inst" "$2/extensions/$n/"
  done
  # the proof lists extension dirs; keep only those
}
run_case() { # name mutate-command
  w=$out/$1/work; s=$out/$1/snap
  copy_tree "$wt" "$w"; copy_tree "$main" "$s"
  (cd "$w" && eval "$2")
  echo "== $1"
  (cd "$w" && "$R" "$wt/dev/conditions-rewrite.R" verify "$s" 2>&1 | \
     grep -E "MISMATCH|files compared")
}
run_case unmutated "true"
# one character of one message in core
run_case mutated_message \
  "sed -i '0,/is not a supported family/s//is not a supportd family/' R/*.R && grep -l 'supportd' R/*.R"
# an argument that is not a message: call. = FALSE dropped at one site
run_case mutated_callarg \
  "sed -i '234s/, call[.] = FALSE)/)/' R/autocor.R && sed -n 234p R/autocor.R"
# a message in an extension
run_case mutated_extension \
  "f=\$(grep -l 'frm_stop(\"' extensions/frmtmb.eam/R/*.R | head -1) && sed -i '0,/frm_stop(\"/s//frm_stop(\"X/' \$f && echo \$f"
rm -rf "$out"

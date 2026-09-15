#!/usr/bin/env bash
# The collision tests on each build, read from the runner outputs.
cd /c/Users/adf44/source/r/frmtmb-wt-samplegen/dev/samplegen-out || exit 1
v() { grep "^$2 " "$1" | awk '{print $2}'; }
row() {
  printf "%-44s %6s %6s %5s %5s %5s\n" "$2" "$(v $1 BLOCKS)" "$(v $1 ASSERT)" \
    "$(v $1 PASS)" "$(v $1 FAIL)" "$(v $1 ERROR)"
  grep "^BAD" "$1" | sed 's/^BAD  */    failing: /'
}
echo '```'
echo "== the collision tests, seen failing, dev/samplegen-runtests*.R =="
printf "%-44s %6s %6s %5s %5s %5s\n" build blocks assert pass fail error
echo "-- extensions/frmtmb.sample/tests/testthat/test-generic-collision.R"
row t-sample-collision-base.txt "BASE core + BASE sample (rellib-r3)"
row t-sample-collision-mutwork.txt "FIX, stancode() generic given work (mutant)"
row t-sample-collision-fix.txt "FIX core + FIX sample"
echo "-- tests/testthat/test-generic-collision.R (core)"
row t-core-collision-base.txt "BASE core (rellib-r3)"
row t-core-collision-fix.txt "FIX core"
echo '```'

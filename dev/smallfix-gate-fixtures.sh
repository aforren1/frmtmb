#!/usr/bin/env bash
# Fixtures for .github/scripts/skip-gate.sh, each with the exit status
# it must give (review finding 12: an empty testthat.Rout passed).
#   bash dev/smallfix-gate-fixtures.sh
# Log: dev/smallfix-out/gate-fixtures.txt
set -u
cd "$(dirname "$0")/.."
F=dev/smallfix-out/gate-fixtures
rm -rf "$F"
sum='[ FAIL 0 | WARN 0 | SKIP 1 | PASS 10 ]'
mk() {
  mkdir -p "$F/$1/x.Rcheck/tests"
  printf '%s' "$2" > "$F/$1/x.Rcheck/tests/testthat.Rout"
}
mk empty ''
mk no-summary '══ Skipped tests (1) ══
• {evd} is not installed (1): '"'"'t.R:1:1'"'"'
'
mk cut-off '> library(testthat)
> test_check("frmtmb")
'
mk clean "$sum
"
mk not-installed "══ Skipped tests (1) ══
• {evd} is not installed (1): 't.R:1:1'
$sum
"
mk cannot-load "══ Skipped tests (1) ══
• {fmesher} cannot be loaded (1): 't.R:1:1'
$sum
"
mkdir -p "$F/none"
fails=0
check() {
  local name=$1 want=$2; shift 2
  out=$(bash .github/scripts/skip-gate.sh "$F/$name" "$@" 2>&1)
  got=$?
  ok=$([ "$got" -eq "$want" ] && echo ok || echo WRONG)
  [ "$ok" = ok ] || fails=$((fails + 1))
  echo "$name (allowed: ${*:-none}): exit $got, want $want: $ok"
  echo "  $(echo "$out" | head -1)"
}
check none 1
check empty 1
check no-summary 1 evd
check cut-off 1
check clean 0
check not-installed 1
check not-installed 0 evd
check cannot-load 1
check cannot-load 1 fmesher
echo "fixtures wrong: $fails"

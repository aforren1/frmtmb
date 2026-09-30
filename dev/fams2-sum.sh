#!/bin/bash
# Sum the RESULT lines of a suite directory, and list any file with a
# failure, error or warning.
f=$1/RESULTS.txt
awk '{for(i=1;i<=NF;i++){split($i,kv,"="); if(kv[1] in s || kv[1] ~ /^(pass|fail|err|skip|warn)$/) s[kv[1]]+=kv[2]}} END{printf "files %d pass %d fail %d err %d skip %d warn %d\n", NR, s["pass"], s["fail"], s["err"], s["skip"], s["warn"]}' "$f"
grep -v "fail=0 err=0 skip=[0-9]* warn=0" "$f"

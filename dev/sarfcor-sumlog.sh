#!/bin/sh
# Per-package totals of a run-par.sh tier log: sh dev/sarfcor-sumlog.sh <log>
awk '/^== .* files ==/{p=$2} /^RESULT/{n[p]++; for(i=3;i<=NF;i++){split($i,a,"="); s[p,a[1]]+=a[2]}} /NO RESULT/{nr[p]++}
END{for(k in n) printf "%s files %d pass %d fail %d err %d skip %d warn %d noresult %d\n", k, n[k], s[k,"pass"], s[k,"fail"], s[k,"err"], s[k,"skip"], s[k,"warn"], nr[k]+0}' "$1" | sort

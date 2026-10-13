#!/usr/bin/env bash
# Lane seflags: summaries pasted into dev/seflags-findings.md.
#   bash dev/seflags-sum.sh scan <run>        per-package totals of a scan
#   bash dev/seflags-sum.sh att <tag>         lane against base, per script
set -u
W=/c/Users/adf44/source/r/frmtmb-wt-seflags
cd "$W"
case "$1" in
  scan)
    s="dev/ciharden-log/$2.sum"
    echo "| package | files | pass | fail | error | skip | warn |"
    echo "|---|---|---|---|---|---|---|"
    awk '$2 == "RESULT" {
      p = $1; n[p]++
      for (i = 4; i <= NF; i++) { split($i, a, "="); t[p, a[1]] += a[2] }
      seen[p] = 1 }
      END { for (p in seen) {
        printf "| %s | %d | %d | %d | %d | %d | %d |\n", p, n[p],
          t[p, "pass"], t[p, "fail"], t[p, "err"], t[p, "skip"], t[p, "warn"]
        N += n[p]; P += t[p, "pass"]; F += t[p, "fail"]; E += t[p, "err"]
        S += t[p, "skip"]; X += t[p, "warn"] }
        printf "| **total** | **%d** | **%d** | **%d** | **%d** | **%d** | **%d** |\n",
          N, P, F, E, S, X }' "$s" | sort
    tail -n 1 "$s"
    grep -c "NO RESULT" "$s" | sed 's/^/files without a RESULT line: /'
    for f in dev/ciharden-log/$2-suite/*.txt; do
      grep -a -h "^lib:" "$f" | head -n 1
    done | awk '{print $NF}' | sort | uniq -c | sed 's/^/frmtmb library: /'
    ;;
  att)
    L="dev/seflags-log/att-$2"
    filt='^lib|^arm|LIB |rellib|seflags-lib|^wrote|BLAS probe|fit s |In addition|no non-missing|In max|^[0-9]+: In|There were|fit [0-9.]+ s'
    for b in "$L"/*-base.txt; do
      s=$(basename "$b" -base.txt)
      [ -f "$L/$s-lane.txt" ] || continue
      d=$(diff <(grep -a -v -E "$filt" "$b") \
               <(grep -a -v -E "$filt" "$L/$s-lane.txt") | grep -c '^[<>]')
      echo "$s: $d differing lines"
    done
    for t in singular smallsd sep; do
      if [ -f "$L/$t-base.tsv" ]; then
        # the verdict columns, not the timings
        # (smallsd column 4 is lme4's own sd, which moves by ulps between
        # processes, dev/rtmb-pitfalls.md item 21)
        k=1-12; [ $t = sep ] && k=1-9; [ $t = smallsd ] && k=1-3,5-11
        n=$(diff <(cut -f$k "$L/$t-base.tsv") <(cut -f$k "$L/$t-lane.tsv") |
              grep -c '^[<>]')
        echo "$t.tsv (columns $k): $n differing lines"
      fi
    done
    ;;
esac

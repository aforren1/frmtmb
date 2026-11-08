#!/usr/bin/env bash
# Lane osa: sum the per-file RESULT lines of the suite runs into
# dev/osa-out/suites.log, per package, with every file that failed,
# errored or warned, every file with a skip, and every file with no
# RESULT line. One awk over all the files, because a process per file is
# slow on this machine under load.
#   bash dev/osa-suites-sum.sh <outdir> [<outdir> ...]
cd /c/Users/adf44/source/r/frmtmb-wt-osa || exit 1
files=()
for d in "$@"; do files+=(dev/osa-out/$d/*--*.txt); done
{
  echo "BEGIN suites"
  awk '
    FNR == 1 {
      if (prev != "" && !hit[prev]) { print "NORESULT", prev; nores++ }
      prev = FILENAME; n = split(FILENAME, a, "/"); base = a[n]
      sub(/--.*/, "", base); pk[FILENAME] = base
    }
    /^lib: frmtmb[ .]/ && !(FILENAME in libl) {
      l = $0; sub(/.*; frmtmb /, "frmtmb ", l); libl[FILENAME] = l; libs[l]++
    }
    /^RESULT/ {
      hit[FILENAME] = 1; p = pk[FILENAME]; files[p]++
      for (i = 2; i <= NF; i++) {
        split($i, kv, "=")
        if (kv[1] == "pass") pass[p] += kv[2]
        if (kv[1] == "fail") fail[p] += kv[2]
        if (kv[1] == "err") err[p] += kv[2]
        if (kv[1] == "skip") skip[p] += kv[2]
        if (kv[1] == "warn") warn[p] += kv[2]
      }
      if ($0 !~ /fail=0 err=0/ || $0 !~ /warn=0/) print "NOTCLEAN", p, $0
      else if ($0 !~ /skip=0 /) print "SKIPS", p, $0
    }
    END {
      if (prev != "" && !hit[prev]) { print "NORESULT", prev; nores++ }
      for (p in files) {
        printf "%-18s files %4d pass %6d fail %d err %d skip %d warn %d\n",
          p, files[p], pass[p], fail[p], err[p], skip[p], warn[p]
      }
      printf "files with no RESULT line: %d\n", nores + 0
      print "frmtmb loaded from (files per lib line):"
      for (l in libs) printf "%6d %s\n", libs[l], l
    }' "${files[@]}"
  echo "END suites"
} > dev/osa-out/suites.log
cat dev/osa-out/suites.log

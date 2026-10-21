#!/usr/bin/env bash
# Release reviewer, 0.73.0: run probes of dev/rel073-rev-probe.R, each
# in its own R process under a 6 GB working-set cap, at most 6 at once.
#   bash dev/rel073-rev-probe.sh <libtag: r11|r10> <probe>...
T=/c/Users/adf44/source/r/frmtmb-wt-release
tag=$1; shift
lib=C:/Users/adf44/source/r/rellib-$tag
for p in "$@"; do
  while [ "$(jobs -rp | wc -l)" -ge 6 ]; do sleep 2; done
  powershell -ExecutionPolicy Bypass -File "$(cygpath -w $T/dev/rel073-rev-capped.ps1)" 6000 \
    "$(cygpath -w $T/dev/rel073-rev-log/probe-$p-$tag.txt)" \
    "$(cygpath -w $T/dev/rel073-rev-probe.R)" "$lib" "$p" &
done
wait
for p in "$@"; do echo "== $p $tag"; cat "$T/dev/rel073-rev-log/probe-$p-$tag.txt"; done

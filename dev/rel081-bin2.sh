#!/usr/bin/env bash
# Bin 2 of the ported brms suite on the final build, rerun by itself:
# chain2's phase 3 ran dev/brmsport-bin2-all.sh in a fresh tree with no
# dev/brmsport-bin2-log/ (gitignored, so the worktree lacks it), and its
# first redirect failed before any R ran. bin 1, the ledger, gen and the
# guards had already run. One R process at a time; it writes only under
# dev/ (R CMD build ignores ^dev$), so it does not touch what a running
# check or session reads. The bin-2 files the pipeline regenerates are
# copied first so a moved hand-written row shows in the diff.
#   bash dev/rel081-bin2.sh > dev/rel081-log/bin2-2.txt 2>&1
set -u
. C:/Users/adf44/source/r/frmtmb-wt-rel081/dev/rel081-env.sh
cd C:/Users/adf44/source/r/frmtmb-wt-rel081 || exit 1
export FRMTMB_PORT_LIB=C:/Users/adf44/source/r/rellib-r20 FRMTMB_PORT_EXTENSIONS=all
export TMP=C:/Users/adf44/AppData/Local/Temp/2 TEMP=C:/Users/adf44/AppData/Local/Temp/2
L=dev/rel081-log
P=$L/bin2-pre
# the copies are the merged tree's, taken once: a rerun compares against
# them, not against its own output
if [ ! -d $P ]; then
  mkdir -p $P
  for f in dev/brmsport-bin2-*.tsv dev/brmsport-bin2-summary.md; do cp "$f" $P/; done
  ( cd dev/brmsport-bin2 && find . -type f | sort | xargs md5sum ) > $P/gen-md5.txt
fi
mkdir -p dev/brmsport-bin2-log
echo "bin2 start $(date +%T)"
sh dev/brmsport-bin2-all.sh; echo "bin2 exit $?"
echo "== regenerated files that moved (CR dropped)"
for f in $P/*.tsv $P/*.md; do
  b=dev/$(basename "$f")
  if ! cmp -s <(tr -d '\r' < "$f") <(tr -d '\r' < "$b"); then
    echo "MOVED $b: $(diff <(tr -d '\r' < "$f") <(tr -d '\r' < "$b") | grep -c '^[<>]') lines"
  fi
done
( cd dev/brmsport-bin2 && find . -type f | sort | xargs md5sum ) > $L/bin2-gen-md5-after.txt
echo "generated bin-2 files moved: $(diff $P/gen-md5.txt $L/bin2-gen-md5-after.txt | grep -c '^>')"
grep -h "^lib:" dev/brmsport-bin2-log/tier-*.txt | sort | uniq -c
echo "BIN2 END $(date +%T)"

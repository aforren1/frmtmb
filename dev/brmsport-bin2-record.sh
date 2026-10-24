#!/bin/sh
# Record every generated bin-2 file, one R process each, into
# dev/brmsport-bin2-log/rec-<pkg>-<topic>.tsv; with "tier" as the first
# argument, run them as the tier would instead (verdicts asserted) into
# tier-<pkg>-<topic>.txt. FRMTMB_PORT_ROOT is the tree (required, as in
# dev/brmsport-record.sh); FRMTMB_PORT_LIB the build (default rellib-r12).
if [ -z "$FRMTMB_PORT_ROOT" ]; then
  echo "FRMTMB_PORT_ROOT is unset: set it to the tree to record."
  exit 1
fi
cd "$FRMTMB_PORT_ROOT" || exit 1
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
mode=record
if [ "${1:-}" = tier ]; then mode=tier; fi
mkdir -p dev/brmsport-bin2-log
for f in dev/brmsport-bin2/*/test-brms-suite-bin2-*.R; do
  pkg=$(basename "$(dirname "$f")")
  topic=$(basename "$f" .R | sed 's/^test-brms-suite-bin2-//')
  if [ "$mode" = record ]; then
    "$R" dev/brmsport-bin2-run.R "$pkg" "$f" \
      "dev/brmsport-bin2-log/rec-$pkg-$topic.tsv" \
      > "dev/brmsport-bin2-log/run-$pkg-$topic.txt" 2>&1
    log="dev/brmsport-bin2-log/run-$pkg-$topic.txt"
  else
    "$R" dev/brmsport-bin2-run.R "$pkg" "$f" \
      > "dev/brmsport-bin2-log/tier-$pkg-$topic.txt" 2>&1
    log="dev/brmsport-bin2-log/tier-$pkg-$topic.txt"
  fi
  grep -h "^lib:" "$log"
  grep -h "^RESULT" "$log" || echo "NO RESULT $pkg $topic"
done

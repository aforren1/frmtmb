#!/usr/bin/env bash
# Is the ported brms suite's record the same at every run? Records every
# generated brms-suite file twice (dev/brmsport-run.R, record mode), one
# R process per file, and lists every assertion whose recorded outcome
# or message differs between the two records.
#
#   bash dev/ciharden-brmsstable.sh <name> <tree> <lib>
#
# <tree>: the checkout whose generated files run (this worktree, or the
# base export dev/ciharden-out/base-tree for the unseeded files).
# <lib>: FRMTMB_PORT_LIB, the library ahead of the user library.
# Output: dev/ciharden-log/brmsstable-<name>/{A,B}/rec-*.tsv and the
# comparison on stdout.
set -u
name="$1"; tree="$2"; lib="$3"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/dev/ciharden-log/brmsstable-$name"
R="/c/Program Files/R/R-4.6.1/bin/Rscript.exe"
export PATH="/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:$PATH"
export FRMTMB_PORT_LIB="$lib"
export FRMTMB_STAN_CACHE="$ROOT/dev/stan-cache"
: "${TMP:?TMP is unset}"
rm -rf "$OUT"; mkdir -p "$OUT/A" "$OUT/B"
cd "$tree"
for arm in A B; do
  for f in tests/testthat/test-brms-suite-*.R \
           extensions/frmtmb.sample/tests/testthat/test-brms-suite-*.R; do
    topic=$(basename "$f" .R | sed 's/^test-brms-suite-//')
    case "$f" in extensions/*) pkg=frmtmb.sample ;; *) pkg=frmtmb ;; esac
    ( "$R" "$ROOT/dev/brmsport-run.R" "$pkg" "$f" \
        "$OUT/$arm/rec-$pkg-$topic.tsv" \
        > "$OUT/$arm/run-$pkg-$topic.txt" 2>&1 ) &
  done
done
wait
cd "$ROOT"
"$R" - "$OUT" <<'EOF'
out <- commandArgs(TRUE)[1]
rd <- function(d) {
  fs <- list.files(file.path(out, d), "^rec-.*[.]tsv$", full.names = TRUE)
  x <- do.call(rbind, lapply(fs, function(f) {
    utils::read.delim(f, header = FALSE, quote = "", colClasses = "character",
                      na.strings = NULL)
  }))
  data.frame(key = paste(x[[1]], x[[2]], x[[3]]), held = x[[5]],
             msg = x[[7]])
}
a <- rd("A"); b <- rd("B")
cat("records:", nrow(a), "and", nrow(b), "rows;",
    length(list.files(file.path(out, "A"), "^rec-")), "and",
    length(list.files(file.path(out, "B"), "^rec-")), "files\n")
m <- merge(a, b, by = "key", all = TRUE, suffixes = c(".a", ".b"))
d <- m[is.na(m$msg.a) | is.na(m$msg.b) | m$held.a != m$held.b |
         m$msg.a != m$msg.b, ]
cat("rows that differ between the two runs:", nrow(d), "\n")
for (i in seq_len(nrow(d))) {
  cat("==", d$key[i], "held", d$held.a[i], "/", d$held.b[i], "\n")
  cat("   A:", substr(d$msg.a[i], 1, 200), "\n")
  cat("   B:", substr(d$msg.b[i], 1, 200), "\n")
}
EOF

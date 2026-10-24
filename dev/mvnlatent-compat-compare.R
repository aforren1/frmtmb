# Compare two dumps of dev/mvnlatent-compat-dump.R as SETS of rows: the
# vocabulary's order moves with the move (contributed features come
# last), so row order is not part of the claim.
#
#   Rscript dev/mvnlatent-compat-compare.R <before.rds> <after.rds> [drop]
#
# With `drop`, rows naming mixture_mvn are left out of BEFORE first:
# the comparison for a session without frmtmb.latent.
a <- commandArgs(trailingOnly = TRUE)
b <- readRDS(a[1])
x <- readRDS(a[2])
key <- function(t) do.call(paste, c(t, sep = "\r"))
# A pair is unordered: its orientation follows the vocabulary order,
# which the move changes for every pair of mixture_mvn with a feature
# that came after it in core's list. Orient each row the same way.
orient <- function(t) {
  if (!all(c("feature_a", "kind_a", "feature_b", "kind_b") %in% names(t))) {
    return(t)
  }
  sw <- t$feature_a > t$feature_b
  fa <- t$feature_a
  ka <- t$kind_a
  t$feature_a[sw] <- t$feature_b[sw]
  t$kind_a[sw] <- t$kind_b[sw]
  t$feature_b[sw] <- fa[sw]
  t$kind_b[sw] <- ka[sw]
  t
}
b$tb <- orient(b$tb)
x$tb <- orient(x$tb)
cat("columns:", paste(names(x$tb), collapse = ", "), "\n")
tb_b <- b$tb
if (length(a) > 2L && a[3] == "drop") {
  hit <- apply(tb_b, 1L, function(r) any(r == "mixture_mvn"))
  cat("before rows naming mixture_mvn, dropped:", sum(hit), "\n")
  tb_b <- tb_b[!hit, , drop = FALSE]
  b$ft <- b$ft[b$ft$name != "mixture_mvn", , drop = FALSE]
}
cat("columns identical:", identical(names(tb_b), names(x$tb)), "\n")
kb <- key(tb_b)
kx <- key(x$tb)
cat("table rows before", length(kb), "after", length(kx), "\n")
cat("rows only before:", sum(!kb %in% kx), "\n")
cat("rows only after:", sum(!kx %in% kb), "\n")
print(head(tb_b[!kb %in% kx, ], 20))
print(head(x$tb[!kx %in% kb, ], 20))
cat("duplicated keys before", sum(duplicated(kb)), "after",
    sum(duplicated(kx)), "\n")
fb <- key(b$ft)
fx <- key(x$ft)
cat("vocabulary before", length(fb), "after", length(fx),
    "only before", sum(!fb %in% fx), "only after", sum(!fx %in% fb), "\n")
print(b$ft[!fb %in% fx, ])
print(x$ft[!fx %in% fb, ])
cat("SETS EQUAL:", setequal(kb, kx) && setequal(fb, fx), "\n")

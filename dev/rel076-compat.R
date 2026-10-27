# The merged compatibility table against 0.75.0, with every extension
# loaded: rellib-r13 (frmtmb 0.75.0 and its eight extensions) against
# rellib-r14 (0.76.0 and its eight). Adapted from dev/rel075-compat.R.
# No lane added a feature or a rule; two lanes reworded hand rules of
# R/compat.R: mocode the mo() row, acprior two group:autocor rows. So
# the pairs and features must be the same, no status may move, and
# every note that moves must have mo() or an autocorrelation term on
# one side. Anything else fails the verdict.
#   Rscript dev/rel076-compat-dump.R dev/rel076-log/compat-base.rds \
#     C:/Users/adf44/source/r/rellib-r13 <the eight>
#   Rscript dev/rel076-compat-dump.R dev/rel076-log/compat-merged.rds \
#     C:/Users/adf44/source/r/rellib-r14 <the eight>
#   Rscript dev/rel076-compat.R > dev/rel076-log/compat.txt
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel076")
a <- readRDS("dev/rel076-log/compat-base.rds")
b <- readRDS("dev/rel076-log/compat-merged.rds")
cat(a$libs, sep = "\n"); cat(b$libs, sep = "\n")
key <- function(t) paste(t$feature_a, t$feature_b, sep = " x ")
ka <- key(a$tab); kb <- key(b$tab)
cat("\npairs before", length(ka), "after", length(kb), "| only before",
    length(setdiff(ka, kb)), "| only after", length(setdiff(kb, ka)), "\n")
fa <- as.data.frame(a$features); fb <- as.data.frame(b$features)
cat("features before", nrow(fa), "after", nrow(fb), "| identical:",
    identical(fa, fb), "\n")
sh <- intersect(ka, kb)
A <- a$tab[match(sh, ka), ]; B <- b$tab[match(sh, kb), ]
st <- names(A)[grepl("status", names(A))][1]
nt <- names(A)[grepl("note", names(A))][1]
ds <- which(A[[st]] != B[[st]])
cat("status differs on a shared pair:", length(ds), "\n")
dn <- which(A[[st]] == B[[st]] & A[[nt]] != B[[nt]])
cat("note differs (status the same):", length(dn), "\n")
is_side <- function(f, k) k == "autocor" | f == "mo()" | grepl("^mo[(]", f)
side <- is_side(A$feature_a[dn], A$kind_a[dn]) |
  is_side(A$feature_b[dn], A$kind_b[dn])
cat("  of them with mo() or an autocorrelation term on one side:",
    sum(side), "\n")
up <- unique(data.frame(before = A[[nt]][dn], after = B[[nt]][dn]))
cat("  distinct note pairs:", nrow(up), "\n")
for (i in seq_len(nrow(up))) {
  cat("  -", substr(up$before[i], 1, 160), "\n    ->",
      substr(up$after[i], 1, 160), "\n")
}
ok <- identical(ka, kb) && identical(fa, fb) && !length(ds) && all(side)
cat("\nVerdict:", if (ok) "PASS" else "FAIL", "\n")

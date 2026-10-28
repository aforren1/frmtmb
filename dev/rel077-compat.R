# The merged compatibility table against 0.76.0, with every extension
# loaded: rellib-r14 (frmtmb 0.76.0 and its eight extensions) against
# rellib-r15 (0.77.0 and its eight). Adapted from dev/rel076-compat.R.
# No lane added a feature or a rule; three lanes reworded hand rules of
# R/compat.R: mointerop the mo() row, sigmainn the residual and se()
# rows of a residual correlation term, reslot the ordinal Intercept row
# of sum-to-zero thresholds. So the pairs and features must be the same,
# no status may move, and every note that moves must have mo(), an
# autocorrelation term or sum_to_zero thresholds on one side. The distinct
# note pairs are listed.
#   Rscript dev/rel077-compat-dump.R dev/rel077-log/compat-base.rds \
#     C:/Users/adf44/source/r/rellib-r14 <the eight>
#   Rscript dev/rel077-compat-dump.R dev/rel077-log/compat-merged.rds \
#     C:/Users/adf44/source/r/rellib-r15 <the eight>
#   Rscript dev/rel077-compat.R > dev/rel077-log/compat.txt
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel077")
a <- readRDS("dev/rel077-log/compat-base.rds")
b <- readRDS("dev/rel077-log/compat-merged.rds")
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
is_side <- function(f, k) k == "autocor" | f == "mo()" | grepl("^mo[(]", f) |
  f == "sum_to_zero"
side <- is_side(A$feature_a[dn], A$kind_a[dn]) |
  is_side(A$feature_b[dn], A$kind_b[dn])
cat("  of them with mo(), an autocorrelation term or sum_to_zero on one side:",
    sum(side), "\n")
up <- unique(data.frame(before = A[[nt]][dn], after = B[[nt]][dn]))
cat("  distinct note pairs:", nrow(up), "\n")
for (i in seq_len(nrow(up))) {
  cat("  -", substr(up$before[i], 1, 160), "\n    ->",
      substr(up$after[i], 1, 160), "\n")
}
ok <- identical(ka, kb) && identical(fa, fb) && !length(ds) && all(side)
cat("\nVerdict:", if (ok) "PASS" else "FAIL", "\n")

# The merged compatibility table against 0.77.0, with every extension
# loaded: rellib-r15 (frmtmb 0.77.0 and its eight extensions) against
# rellib-r16 (0.78.0 and its eight). Adapted from dev/rel077-compat.R.
# No lane added a feature or a rule; two lanes reworded hand rules: meme
# the me() x prior row of R/compat.R (conditional to works, brms's
# me() prior classes) and rescorpp frmtmb.sample's rescor x frm_sample
# note (R/zzz.R: the joint draw and density). So the pairs and features
# must be the same; a status may move only on a pair with me() on one
# side and prior on the other, and a note only on such a pair or on
# rescor x frm_sample. Every pair that moves is listed.
#   Rscript dev/rel078-compat-dump.R dev/rel078-log/compat-base.rds \
#     C:/Users/adf44/source/r/rellib-r15 <the eight>
#   Rscript dev/rel078-compat-dump.R dev/rel078-log/compat-merged.rds \
#     C:/Users/adf44/source/r/rellib-r16 <the eight>
#   Rscript dev/rel078-compat.R > dev/rel078-log/compat.txt
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel078")
a <- readRDS("dev/rel078-log/compat-base.rds")
b <- readRDS("dev/rel078-log/compat-merged.rds")
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
pair_is <- function(i, x, y) {
  (A$feature_a[i] == x & A$feature_b[i] == y) |
    (A$feature_a[i] == y & A$feature_b[i] == x)
}
expected <- function(i) pair_is(i, "me()", "prior") |
  pair_is(i, "rescor", "frm_sample")
ds <- which(A[[st]] != B[[st]])
cat("status differs on a shared pair:", length(ds), "\n")
for (i in ds) cat("  -", sh[i], ":", A[[st]][i], "->", B[[st]][i], "\n")
dn <- which(A[[st]] == B[[st]] & A[[nt]] != B[[nt]])
cat("note differs (status the same):", length(dn), "\n")
for (i in dn) cat("  -", sh[i], "\n")
up <- unique(data.frame(before = A[[nt]][c(ds, dn)],
                        after = B[[nt]][c(ds, dn)]))
cat("distinct note pairs:", nrow(up), "\n")
for (i in seq_len(nrow(up))) {
  cat("  -", substr(up$before[i], 1, 200), "\n    ->",
      substr(up$after[i], 1, 200), "\n")
}
ok <- identical(ka, kb) && identical(fa, fb) &&
  all(expected(ds)) && all(expected(dn)) &&
  all(B[[st]][ds] == "works")
cat("\nVerdict:", if (ok) "PASS" else "FAIL", "\n")

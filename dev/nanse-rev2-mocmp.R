# Reviewer: compare the trial merge's mo() TSV with the lane's p1b TSV
# and round 0's, column by column.
#   Rscript dev/nanse-rev2-mocmp.R
w <- "C:/Users/adf44/source/r/frmtmb-wt-nanse/dev/"
a <- read.delim(paste0(w, "nanse-log/mo-lane-p1b.tsv"))
b <- read.delim(paste0(w, "nanse-rev2-log/mo-merge.tsv"))
z <- read.delim(paste0(w, "nanse-log/mo-lane-final.tsv"))
stopifnot(identical(a$seed, b$seed), identical(a$form, b$form))
for (cn in intersect(names(a), names(b))) {
  d <- which(!(a[[cn]] == b[[cn]] | (is.na(a[[cn]]) & is.na(b[[cn]]))) |
               xor(is.na(a[[cn]]), is.na(b[[cn]])))
  if (length(d)) cat("p1b vs merge:", cn, "differs on", length(d), "rows; seeds",
                     paste(head(paste(a$form[d], a$seed[d]), 8), collapse = " "), "\n")
}
cat("ce_ok false: round0", sum(!z$ce_ok), " p1b", sum(!a$ce_ok), " merge",
    sum(!b$ce_ok), "\n")
d <- which(a$ce_ok != z$ce_ok)
cat("round0 vs p1b ce_ok differ on:", paste(a$form[d], a$seed[d]), "\n")

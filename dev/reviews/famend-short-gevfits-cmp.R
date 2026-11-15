# Re-check of lane famend after fix round 2: compare the GEV fits of
# famend-short-gevfits.R on the base and the lane builds.
#   Rscript dev/reviews/famend-short-gevfits-cmp.R <base.tsv> <lane.tsv>
a <- commandArgs(trailingOnly = TRUE)
b <- utils::read.delim(a[1], colClasses = "character")
l <- utils::read.delim(a[2], colClasses = "character")
stopifnot(identical(b[, 1:3], l[, 1:3]))
ob <- as.numeric(b$objective)
ol <- as.numeric(l$objective)
out <- as.numeric(b$n_outside) > 0 | as.numeric(l$n_outside) > 0
est_b <- lapply(strsplit(b$est, ","), as.numeric)
est_l <- lapply(strsplit(l$est, ","), as.numeric)
maxrel <- mapply(function(u, v) {
  if (length(u) != length(v)) return(Inf)
  max(abs(u - v) / pmax(abs(u), 1e-8))
}, est_b, est_l)
same_obj <- ob == ol
cat(sprintf("fits %d; with an outside trial point (either build) %d\n",
            nrow(b), sum(out)))
cat(sprintf("code equal %d of %d\n", sum(b$code == l$code), nrow(b)))
cat(sprintf("objective identical %d; lane lower %d; lane higher %d; NA %d\n",
            sum(same_obj, na.rm = TRUE), sum(ol < ob, na.rm = TRUE),
            sum(ol > ob, na.rm = TRUE), sum(is.na(ob) | is.na(ol))))
cat(sprintf("estimates identical to 8 digits %d; max relative change %.3g\n",
            sum(b$est == l$est), max(maxrel[is.finite(maxrel)])))
d <- which(!(same_obj %in% TRUE) | b$est != l$est | b$code != l$code)
if (length(d)) {
  cat("rows that differ:\n")
  print(data.frame(xi0 = b$xi0[d], seed = b$seed[d], form = b$form[d],
                   code_b = b$code[d], code_l = l$code[d],
                   obj_b = b$objective[d], obj_l = l$objective[d],
                   d_obj = ol[d] - ob[d], est_rel = signif(maxrel[d], 3),
                   out_b = b$n_outside[d], out_l = l$n_outside[d]),
        row.names = FALSE)
}

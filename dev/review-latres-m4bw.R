# The two unmeasured families, from the benign set.
x <- do.call(rbind, lapply(list.files("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/m4b2", "tsv$", full.names = TRUE),
                           function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
x <- x[is.finite(x$ll), ]
for (d in c("bb_post_lo", "zoib_post_lo", "bb_post_b", "zoib_post_b")) {
  s <- x[x$design == d, ]
  ok <- abs(s$ls - s$x_ls) < log(1.25) & abs(s$ar - s$x_ar) < .3 & abs(s$ll - s$x_at_frm) < 1
  cat(d, "fits", nrow(s), "code", sum(s$code != 0), "nwarn", sum(s$nwarn > 0), "sderr", round(median(exp(s$ls)), 3), round(median(exp(s$x_ls)), 3),
      "ar", round(median(s$ar), 2), round(median(s$x_ar), 2), "gap", round(median(s$ll - s$x_at_frm), 2), "true", sum(ok), "\n")
}

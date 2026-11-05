# Is skew_normal()'s coverage of 0.73 about Laplace? The exact MLE's own distance
# from the truth, in the fit's SE, and the correction there.
for (d in c("sn_u1", "sn_u10", "binom4", "nb_m1", "nb_s1")) {
  x <- read.table(paste0("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/tk/m1b-", d, ".tsv"), sep = "\t", header = TRUE, quote = "\"", comment.char = "")
  t0 <- log(x$s_true[1])
  se <- x$se_ls
  cat(sprintf("%-8s fits %d  cover frm %.3f  cover exact MLE with frm's SE %.3f  median |ls - exact| %.3f  median |C| %.3f  median |shift| %.3f  wrong per fit %.2f  SE NaN %d\n", d, nrow(x),
              mean(abs(x$ls - t0) <= 1.96 * se, na.rm = TRUE), mean(abs(x$x_ls - t0) <= 1.96 * se, na.rm = TRUE), median(abs(x$ls - x$x_ls), na.rm = TRUE), median(abs(x$C), na.rm = TRUE), median(pmax(abs(x$b_ls), abs(x$b_ar)), na.rm = TRUE),
              mean(abs(x$ls - x$x_ls) > log(1.25) | abs(x$ar - x$x_ar) > .3, na.rm = TRUE), sum(!is.finite(se))))
}

# The gate on the negbinomial designs, against wrong per fit and the lane's warning.
fs <- list.files("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/tk2", "(nb|nb_lo|nb_m1|nb_m5|nb_m20|nb_s1|f_nb1|f_nb8)[.]tsv$", full.names = TRUE)
for (f in fs) {
  x <- read.table(f, sep = "\t", header = TRUE, quote = "\"", comment.char = "")
  s <- pmax(abs(x$b_ls), abs(x$b_ar)); s[!is.finite(s)] <- Inf
  w <- abs(x$ls - x$x_ls) > log(1.25) | abs(x$ar - x$x_ar) > .3
  cat(sprintf("%-8s fits %d wrong %d lane %d gate %d\n", x$design[1], nrow(x), sum(w, na.rm = TRUE), sum(x$warned == 1), sum(s > .15)))
}

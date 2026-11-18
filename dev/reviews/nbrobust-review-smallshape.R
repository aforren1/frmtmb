# Review of lane nbrobust: value error at a small shape, nb_lpmf() against
# dnbinom_robust() and stats::dnbinom(), from
# dev/reviews/nbrobust-review-formula.tsv
t <- read.delim("dev/reviews/nbrobust-review-formula.tsv")
for (ls in c(-10, -3, 0)) {
  a <- t[t$arm == "scalar" & t$lsh == ls, ]
  b <- t[t$arm == "robust" & t$lsh == ls, ]
  stopifnot(identical(a$y, b$y), identical(a$lmu, b$lmu))
  out <- data.frame(y = a$y, lmu = a$lmu, vref = signif(a$vref, 6),
                    nb_abs = signif(a$absv, 3), robust_abs = signif(b$absv, 3),
                    stats_abs = signif(abs(a$vstats - a$vref), 3),
                    nb_units = signif(a$absv / (.Machine$double.eps *
                      (abs(a$vref) + a$y * abs(ls) + 1)), 3))
  cat("log shape", ls, ": rows where nb_lpmf's error exceeds 4x the larger",
      "of the other two and 1e-13\n")
  w <- out$nb_abs > 4 * pmax(out$robust_abs, out$stats_abs) &
    out$nb_abs > 1e-13
  print(out[w, ], row.names = FALSE)
  cat("max abs: nb", max(out$nb_abs), " robust", max(out$robust_abs),
      " stats", max(out$stats_abs), "\n\n")
}

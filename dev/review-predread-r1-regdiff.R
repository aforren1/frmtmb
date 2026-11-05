# Short check, item 3 of the first review: which columns moved on a model with no me().
a <- readRDS("dev/review-predread-out/regress-before.rds")
b <- readRDS("dev/review-predread-out/regress-after.rds")
for (k in c("poisson|fitted", "poisson|fitted_nd", "zip|fitted", "bernoulli|fitted_nd", "lognormal|fitted_nd", "hurdle|fitted_nd")) {
  x <- a[[k]]
  y <- b[[k]]
  cat("\n", k, "\n")
  print(signif(rbind(before = x[1, ], after = y[1, ]), 5))
  cat("Estimate identical:", identical(x[, "Estimate"], y[, "Estimate"]), " Est.Error identical:", identical(x[, "Est.Error"], y[, "Est.Error"]), "\n")
  cat("Q columns, largest rel:", signif(max(abs(x[, c("Q2.5", "Q97.5")] - y[, c("Q2.5", "Q97.5")]) / abs(x[, c("Q2.5", "Q97.5")])), 3), "\n")
}
for (k in c("lognormal|CE", "hurdle|CE", "ar_cond|CE_predict")) {
  x <- a[[k]]
  y <- b[[k]]
  cat("\n", k, " panels", names(x), "\n")
  print(signif(cbind(as.data.frame(x[[1L]]), as.data.frame(y[[1L]])), 4))
}

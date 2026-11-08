source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
ns <- asNamespace("frmtmb")
set.seed(501); n <- 200
x1 <- rnorm(n); invisible(rnorm(n)); invisible(rnorm(n)); invisible(rnorm(n))
xx <- runif(n, 1000, 1010)
for (deg in 5:8) {
  X <- model.matrix(~ poly(xx, deg, raw = TRUE))
  s <- sqrt(colSums(X^2) / n)
  q <- qr(X / rep(s, each = n), tol = 1e-11)
  R <- qr.R(q)
  Tb <- tryCatch(ns$design_whitener(X, 5e6), error = function(e) e)
  cat("deg", deg, "frame rank@1e-11", qr(X, tol = 1e-11)$rank,
      "scaled rank@1e-11", q$rank, "rcond kappa", format(1 / rcond(R, triangular = TRUE), digits = 3),
      "whitener", if (is.null(Tb)) "NULL" else if (inherits(Tb, "error")) conditionMessage(Tb) else
        paste("ok finite", all(is.finite(Tb)), "inv finite", all(is.finite(attr(Tb, "inv")))), "\n")
  cat("   |diag R| min/max", format(min(abs(diag(R))) / max(abs(diag(R))), digits = 3), "\n")
}

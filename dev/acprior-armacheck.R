# frmtmb's cov = TRUE ar(1), ma(1) and arma(1, 1) correlation against
# brms's cholesky_cor_*() covariance divided by its diagonal, at the same
# natural coefficients; and the diagonal itself, the ratio of brms's
# marginal variance to its sigma^2. Log: dev/acprior-out/armacheck.txt
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
d <- data.frame(t = rep(1:5, 10), g = rep(1:10, each = 5), y = rnorm(50))
brms_cov <- function(kind, ar, ma, n) {
  if (kind == "ar") {
    m <- ar^abs(outer(1:n, 1:n, "-")) / (1 - ar^2)
  } else if (kind == "ma") {
    m <- diag(1 + ma^2, n)
    m[abs(row(m) - col(m)) == 1] <- ma
  } else {
    g1 <- (1 + ar * ma) * (ar + ma)
    m <- diag(1 + ma^2 + 2 * ar * ma, n)
    for (i in 2:n) for (j in 1:(i - 1)) {
      m[i, j] <- m[j, i] <- g1 * ar^(i - j - 1)
    }
    m <- m / (1 - ar^2)
  }
  m
}
for (kind in c("ar", "ma", "arma")) {
  fo <- switch(kind, ar = y ~ ar(t, g, cov = TRUE),
               ma = y ~ ma(t, g, cov = TRUE),
               arma = y ~ arma(t, g, cov = TRUE))
  fr <- frm(fo, data = d, dry_run = "frame")
  ac <- fr$autocor[[1]]
  for (rep in 1:3) {
    th <- rnorm(ac$npar)
    nat <- as.numeric(frmtmb:::autocor_natural(th, ac))
    ar <- if (kind != "ma") nat[1] else 0
    ma <- if (kind == "ma") nat[1] else if (kind == "arma") nat[2] else 0
    R <- as.matrix(frmtmb:::autocor_cor(th, ac))
    B <- brms_cov(kind, ar, ma, ac$d)
    cat(sprintf("%-4s ar %+.4f ma %+.4f  max|R - cov2cor(B)| %.2e  diag(B) %.6f\n",
                kind, ar, ma, max(abs(R - stats::cov2cor(B))), B[1, 1]))
  }
}

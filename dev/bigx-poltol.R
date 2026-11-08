# Lane bigx, fix round, major 2: which alias tolerance on the scaled
# design leaves only columns a fit can estimate. For each tolerance,
# raw polynomials of x on [1000, 1010] (seeds 501 to 503, degree 5 to 8)
# and a shifted column x + m (poisson with (1 | g), m = 1e8 and 1e9, data
# seeds 1 and 2): columns kept, the kept design's scaled condition
# number, the code, finite SEs, and the log-likelihood against glm() on
# the same kept columns (or the centered fit for the shifted column).
#   BIGX_ARM=after Rscript dev/bigx-poltol.R
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
ns <- asNamespace("frmtmb")
for (tol in c(1e-11, 1e-10, 1e-9)) {
  assignInNamespace("design_alias_tol", tol, "frmtmb")
  cat("== tol", tol, "\n")
  for (s in 501:503) for (deg in 5:8) {
    set.seed(s)
    n <- 200
    x1 <- rnorm(n); invisible(rnorm(n)); invisible(rnorm(n)); invisible(rnorm(n))
    xx <- runif(n, 1000, 1010)
    d <- data.frame(xx = xx, yy = sin(xx / 2) + rnorm(n, 0, 0.3))
    fo <- as.formula(sprintf("yy ~ poly(xx, %d, raw = TRUE)", deg))
    x <- suppressMessages(bigx_fit(bf(fo), family = gaussian(), data = d))
    if (x$err) { cat("  ERROR", conditionMessage(x$f), "\n"); next }
    f <- x$f
    X <- f$frame$linpreds[[1]]$X
    kap <- ns$lp_design_kappa(list(X = X))
    g <- glm.fit(X, d$yy)
    llg <- -n / 2 * (log(2 * pi * sum(g$residuals^2) / n) + 1)
    se <- tryCatch(suppressWarnings(sqrt(diag(vcov(f)))), error = function(e) NA)
    cat(sprintf("  seed %d deg %d kept %d kappa %.2e code %d ll - glm(kept) %+.2e SEs finite %d/%d warn %d\n",
                s, deg, ncol(X), kap, f$opt$convergence,
                as.numeric(logLik(f)) - llg, sum(is.finite(se)), length(se),
                length(x$w)))
  }
  for (m in c(1e8, 1e9)) for (s in 1:2) {
    d <- bigx_data("poisson", m, s)
    fc <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
    x <- suppressMessages(bigx_fit(bf(y ~ w + (1 | g)), family = poisson(),
                                   data = d))
    f <- x$f
    cat(sprintf("  shift %g seed %d kept %d code %d dll %+.2e warn %d\n", m, s,
                length(f$estimates[["beta"]]), f$opt$convergence,
                as.numeric(logLik(f)) - as.numeric(logLik(fc)), length(x$w)))
  }
}

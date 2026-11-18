# Review of lane sparsese: how often the whole-inverse sparse route
# (jcs_joint_cov(white = TRUE, whole = TRUE)) returns NULL, so that
# joint_cov_whole() falls back to the natural inverse, on lost-sd fits:
# the first test of test-se-reparam-reml.R's data (200 rows, 20 groups,
# degree 7, no group effect), 40 ulp perturbations, limit 20.
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
q <- function(expr) suppressWarnings(suppressMessages(expr))
options(frmtmb.joint_root_max = 20)
res <- NULL
for (seed in 1:4) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  y0 <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  X <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, d)
  for (k in 0:9) {
    d$y <- y0 * (1 + 2e-15 * k)
    fit <- q(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d,
                 REML = TRUE))
    sdr <- frmtmb:::sdr_of(fit)
    sp <- frmtmb:::jcs_joint_cov(fit, white = TRUE, whole = TRUE)
    V <- q(vcov(fit))
    s2 <- exp(2 * fit$opt$par[["betad"]])
    # with the sd lost the GLS weights are 1: s2 (X'X)^-1, via QR
    qx <- qr(X / rep(sqrt(colSums(X^2)), each = nrow(X)))
    Ri <- matrix(0, ncol(X), ncol(X))
    Ri[qx$pivot, ] <- backsolve(qr.R(qx), diag(ncol(X))) /
      sqrt(colSums(X^2))[qx$pivot]
    ols <- sqrt(s2 * rowSums(Ri^2))
    rr <- range(sqrt(diag(V)) / ols - 1)
    res <- rbind(res, data.frame(seed = seed, k = k,
                                 lost = paste(sdr$se_lost, collapse = ""),
                                 theta = fit$opt$par[["theta"]],
                                 whole_null = is.null(sp),
                                 vcov_lo = rr[1], vcov_hi = rr[2]))
  }
}
print(res, digits = 3)
cat("lost fits:", sum(res$lost != ""), " whole route NULL:",
    sum(res$whole_null), "\n")

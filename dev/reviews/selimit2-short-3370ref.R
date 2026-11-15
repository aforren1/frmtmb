# Short review of lane selimit2, F4: as selimit2-review-3370ref.R, any seed.
# The conditional SE of the cost fit of selimit2-review-3370.R (`ng` groups
# of 2 rows, degree 7) by GLS at the fitted group variance in an orthogonal
# basis; a fitted SE (which adds the variance parameters' propagation) must
# not fall below it.
#   Rscript selimit2-short-3370ref.R <base|lane> <groups> <seed>
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
ng <- as.integer(a[2])
seed <- as.integer(a[3])
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
set.seed(seed)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                f = factor(rep(c("a", "b"), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fit <- withCallingHandlers(suppressMessages(frm(
  y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d, REML = TRUE)),
  warning = function(w) {
    cat("fit W:", substr(conditionMessage(w), 1, 200), "
")
    invokeRestart("muffleWarning")
  })
est <- fit$opt$par
cat("par:", paste(names(est), format(est, digits = 6)), "\n")
lost <- frmtmb:::sdr_of(fit)$se_lost
cat("lost:", paste(names(lost), lost), "\n")
s2 <- exp(2 * fit$estimates$betad[1])
t <- exp(2 * (fit$estimates$theta[1] - fit$estimates$betad[1]))
pb <- stats::poly(d$x, 7)
X <- cbind(1, pb)
nd <- data.frame(x = c(1, 1.5, 2))
A <- cbind(1, stats::predict(pb, nd$x))
# whitened GLS: subtract a per group, as test-se-reparam-reml.R's reml_gls
Xt <- X
for (ix in split(seq_len(n), d$g)) {
  m <- length(ix)
  aa <- (1 - 1 / sqrt(1 + m * t)) / m
  Xt[ix, ] <- X[ix, , drop = FALSE] - aa * rep(colSums(X[ix, , drop = FALSE]),
                                               each = m)
}
Ri <- backsolve(qr.R(qr(Xt)), diag(ncol(Xt)))
cat("variance ratio t", format(t, digits = 4), "\n")
cat("conditional SE:", format(sqrt(s2 * rowSums((A %*% Ri)^2)), digits = 7),
    "\n")
cat("OLS SE (t = 0):", format(sqrt(s2 * rowSums((A %*% backsolve(
  qr.R(qr(X)), diag(ncol(X))))^2)), digits = 7), "\n")
nd2 <- data.frame(x = c(1, 1.5, 2), g = d$g[1])
for (k in 1:7) nd2[[paste0("x", k)]] <- nd2$x^k
fp <- withCallingHandlers(stats::fitted(fit, newdata = nd2, re_formula = NA),
  warning = function(w) {
    cat("fitted W:", substr(conditionMessage(w), 1, 200), "
")
    invokeRestart("muffleWarning")
  })
cat("fitted SE:", format(fp[, "Est.Error"], digits = 7), "\n")

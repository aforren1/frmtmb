# Reviewer of lane seflags, final check of the second defect: degree-7
# raw polynomial with (1 | g), seed 7 (nothing lost, a span
# re-measured). fitted(newdata, re_formula = NA) at x = 1, 1.5, 2, 2.5
# against the exact prediction SE x0' (X' V^-1 X)^-1 x0 at the fit's own
# sigma and sd, from a QR of the whitened design; whether the lane kept
# a repaired covariance.
#   Rscript dev/seflags-rev4-deg7pred.R <lane|base|op3|m9>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-common.R")
set.seed(7)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
  stats::rnorm(200, 0, 0.2)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:7, collapse = " + ")
r <- cap(frm(stats::as.formula(paste("y ~", rhs, "+ (1 | g)")), data = d))
f <- r$v
sdr <- suppressWarnings(ns$sdr_of(f))
on <- ns$outer_par_names(f)
nd <- data.frame(x = c(1, 1.5, 2, 2.5), g = d$g[1])
for (k in 1:7) nd[[paste0("x", k)]] <- nd$x^k
p <- cap(fitted(f, newdata = nd, re_formula = NA))
X <- stats::model.matrix(stats::as.formula(paste("~", rhs)), d)
X0 <- stats::model.matrix(stats::as.formula(paste("~", rhs)), nd)
Z <- stats::model.matrix(~ 0 + g, d)
s <- exp(f$opt$par[grepl("^sigma", on)])
t <- exp(f$opt$par[on == "theta_1"])
L <- chol(s^2 * diag(nrow(X)) + t^2 * Z %*% t(Z))
q <- qr(backsolve(L, X, transpose = TRUE), tol = 1e-30)
A <- backsolve(qr.R(q), t(X0[, q$pivot, drop = FALSE]), transpose = TRUE)
ex <- sqrt(colSums(A^2))
cat(sprintf("code %d | lost %s | fit warnings %d | repaired covariance %s\n",
            f$opt$convergence, lost_of(f), length(r$w),
            !is.null(sdr$cov_fixed_prop)))
cat("fitted se ", format(p$v[, 2], digits = 5), "| warnings", length(p$w), "\n")
cat("exact     ", format(ex, digits = 5), "\n")
cat("ratio     ", format(p$v[, 2] / ex, digits = 5), "\n")

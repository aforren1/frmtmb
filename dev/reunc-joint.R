# Lane wt-reunc: is predict(summary = FALSE) JOINTLY right across rows?
#
# Gaussian one-way random intercept, so the conditional law of b given
# the data and the parameters is analytic:
#   Var(b_g | y, beta, sigma, tau) = sigma^2 tau^2 / (sigma^2 + n_g tau^2)
# and two rows of one group, drawn at the estimates, correlate at
#   v_g / (v_g + sigma^2),
# while rows of two different groups do not correlate at all. With the
# parameter draw on, the model-implied covariance of two rows is
# a_i' V a_j from frm_lp_basis(), the joint covariance of (beta, b).
#
#   Rscript dev/reunc-joint.R <lib> [mode]
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
mode <- if (length(a) > 1) a[2] else "ML"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), " mode", mode, "\n")
set.seed(4)
G <- 6; m <- 4
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 1)[d$g] + rnorm(G * m, 0, 0.8)
args <- list(bf(y ~ x + (1 | g)), data = d)
if (mode == "REML") args$REML <- TRUE
if (mode == "profile") args$control <- frmtmb_control(profile = TRUE)
fit <- do.call(frm, args)
s2 <- sigma(fit)^2
t2 <- exp(2 * fit$estimates[["theta"]][1])
v <- s2 * t2 / (s2 + m * t2)
nd <- data.frame(x = c(0, 1, 0, -1), g = factor(c(1, 1, 2, 3),
                                                levels = levels(d$g)))
N <- 20000

# The 2026-09-23 rename made `propagate_error = FALSE` hold the group
# effects as well as the parameters, so there is no public call left
# that draws b alone. This rebuilds exactly what the old flag gave:
# the drawer at the estimates (delta = NULL is "no parameter
# deviation"), pushed through the linear predictor, plus the noise.
plug_draws <- function(fit, nd, N, seed) {
  set.seed(seed)
  drawer <- frmtmb:::predict_b_drawer(fit, NULL, NULL, N, TRUE)
  stopifnot(!is.null(drawer))
  s2 <- sigma(fit)^2
  t(vapply(seq_len(N), function(s) {
    fs <- drawer(fit, s)
    as.numeric(frmtmb:::fitted_point(fs, nd)) +
      stats::rnorm(nrow(nd), 0, sqrt(s2))
  }, numeric(nrow(nd))))
}
P <- plug_draws(fit, nd, N, 1)
C <- cov(P)
want_same <- v
cat(sprintf("plug-in, same group:   cov %.5f  analytic %.5f  ratio %.4f\n",
            C[1, 2], want_same, C[1, 2] / want_same))
cat(sprintf("plug-in, var row 1:    %.5f  analytic %.5f  ratio %.4f\n",
            C[1, 1], v + s2, C[1, 1] / (v + s2)))
cat(sprintf("plug-in, other group:  cov %.5f %.5f  (analytic 0; mc se %.5f)\n",
            C[1, 3], C[3, 4], sqrt(C[1, 1] * C[3, 3] / N)))
set.seed(2)
P2 <- predict(fit, newdata = nd, summary = FALSE, ndraws = N)
C2 <- cov(P2)
lb <- frm_lp_basis(fit, newdata = nd)
M <- lb$A %*% lb$V %*% t(lb$A)
for (ij in list(c(1, 2), c(1, 3), c(3, 4))) {
  i <- ij[1]; j <- ij[2]
  se <- sqrt((C2[i, i] * C2[j, j] + C2[i, j]^2) / N)
  cat(sprintf("full, cov[%d,%d]: %.5f  implied a_i'Va_j %.5f  z %.2f\n",
              i, j, C2[i, j], M[i, j], (C2[i, j] - M[i, j]) / se))
}
for (i in 1:4) {
  cat(sprintf("full, var[%d]: %.5f  implied a'Va + sigma^2 %.5f\n", i,
              C2[i, i], M[i, i] + s2))
}
fl <- fitted(fit, newdata = nd)
cat("fitted Est.Error^2:", round(fl[, 2]^2, 5), "\n")
cat("a'Va              :", round(diag(M), 5), "\n")

# Lane wt-reunc: is the plug-in variance of one row exactly
# Var(b_g | y) + sigma^2, or is the 0.979 ratio at 20000 draws
# (dev/reunc-log/joint-lane.txt) something other than Monte Carlo
# error? Three seeds at 100000 draws each, on the same fit.
# The library is an ARGUMENT with a default, not a constant: pinned to
# one session's private library this script silently ran against a
# stale build after a later session rebuilt the package elsewhere.
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
set.seed(4)
G <- 6; m <- 4
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 1)[d$g] + rnorm(G * m, 0, 0.8)
fit <- frm(bf(y ~ x + (1 | g)), data = d)
s2 <- sigma(fit)^2
t2 <- exp(2 * fit$estimates[["theta"]][1])
v <- s2 * t2 / (s2 + m * t2)
nd <- data.frame(x = c(0, 1), g = factor(c(1, 1), levels = levels(d$g)))
N <- 100000

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
for (sd in 1:3) {
  P <- plug_draws(fit, nd, N, sd)
  C <- cov(P)
  cat(sprintf(paste0("seed %d: var ratio %.4f (z %.2f), same-group cov ",
                     "ratio %.4f (z %.2f)\n"), sd, C[1, 1] / (v + s2),
              (C[1, 1] - (v + s2)) / (C[1, 1] * sqrt(2 / N)), C[1, 2] / v,
              (C[1, 2] - v) / sqrt((C[1, 1] * C[2, 2] + C[1, 2]^2) / N)))
}

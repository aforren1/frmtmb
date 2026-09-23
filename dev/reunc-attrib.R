# Lane wt-reunc: what is left after predict() carries the group
# effects' conditional law.
#
# dev/reunc-coverage.R's oracle arm builds the exact interval at the
# TRUE sigma and tau and covers 0.95 by construction. predict() has to
# use the fit's ESTIMATES, and with few groups tau is estimated badly.
# This arm is the oracle with the fit's own estimates plugged in, and
# nothing else changed: no draws, no frmtmb prediction code. If its
# coverage lands where predict()'s lands, the remaining shortfall is
# the variance parameters and not the construction.
#
# Same designs, same seeds and same data as dev/reunc-coverage.R.
#
#   Rscript dev/reunc-attrib.R <lib> <design> <nrep>
a <- commandArgs(trailingOnly = TRUE)
lib <- a[1]
design <- a[2]
nrep <- as.integer(a[3])
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressPackageStartupMessages(library(frmtmb))
ng <- switch(design, mixed12 = 12L, few4 = 4L)
n <- switch(design, mixed12 = 60L, few4 = 40L)
seed0 <- switch(design, mixed12 = 20260922L, few4 = 20260923L)
m <- 10L
TAU <- 0.7
SIGMA <- 1
BETA <- c(1, 0.5)
z <- stats::qnorm(0.975)

hend <- function(d, nd, s2, t2) {
  X <- cbind(1, d$x)
  Z <- stats::model.matrix(~ 0 + g, d)
  C <- rbind(cbind(crossprod(X), crossprod(X, Z)),
             cbind(crossprod(Z, X), crossprod(Z) + diag(s2 / t2, ncol(Z)))) /
    s2
  Ci <- solve(C)
  est <- drop(Ci %*% (c(crossprod(X, d$y), crossprod(Z, d$y)) / s2))
  An <- cbind(1, nd$x, stats::model.matrix(~ 0 + g, nd))
  list(fit = drop(An %*% est), pev = rowSums((An %*% Ci) * An))
}
hit <- function(h, y, s2) {
  w <- z * sqrt(s2 + h$pev)
  mean(y >= h$fit - w & y <= h$fit + w)
}
arms <- c("true", "ml_est", "reml_est")
acc <- stats::setNames(rep(list(numeric(0)), length(arms)), arms)
for (r in seq_len(nrep)) {
  set.seed(seed0 + 1000L * r)
  d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng),
                                               length.out = n)))
  u <- rnorm(ng, 0, TAU)
  d$y <- rnorm(n, BETA[1] + BETA[2] * d$x + u[d$g], SIGMA)
  gsel <- sample(seq_len(ng), m, TRUE)
  nd <- data.frame(x = rnorm(m), g = factor(gsel, levels = levels(d$g)))
  y_new <- rnorm(m, BETA[1] + BETA[2] * nd$x + u[gsel], SIGMA)
  f_ml <- tryCatch(suppressWarnings(frm(bf(y ~ x + (1 | g)), data = d)),
                   error = function(e) NULL)
  f_re <- tryCatch(suppressWarnings(frm(bf(y ~ x + (1 | g)), data = d,
                                        REML = TRUE)),
                   error = function(e) NULL)
  if (is.null(f_ml) || is.null(f_re)) next
  s2ml <- sigma(f_ml)^2
  t2ml <- exp(2 * f_ml$estimates[["theta"]][1])
  s2re <- sigma(f_re)^2
  t2re <- exp(2 * f_re$estimates[["theta"]][1])
  acc$true <- c(acc$true, hit(hend(d, nd, SIGMA^2, TAU^2), y_new, SIGMA^2))
  acc$ml_est <- c(acc$ml_est, hit(hend(d, nd, s2ml, t2ml), y_new, s2ml))
  acc$reml_est <- c(acc$reml_est, hit(hend(d, nd, s2re, t2re), y_new, s2re))
}
cat(sprintf("design %s, %d replicates of %d points, script dev/reunc-attrib.R\n",
            design, length(acc$true), m))
for (nm in arms) {
  v <- acc[[nm]]
  se <- stats::sd(v) / sqrt(length(v))
  cat(sprintf("%-9s coverage %.4f  rep se %.4f  z %.2f\n", nm, mean(v), se,
              (mean(v) - 0.95) / se))
}
# and how far the variance parameters are from the truth
cat("\n")

# REML at low noise against its closed form. Seed 20261102.
#   Rscript dev/ncgp-reml.R <lib or base>
a <- commandArgs(TRUE)
libs <- c("C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (a[1] != "base") libs <- c(a[1], libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
for (s in c(1e-1, 1e-2, 1e-3)) {
  set.seed(20261102)
  n0 <- 60L
  x <- round(runif(n0, 0, 6), 3)
  e0 <- rnorm(n0)
  d <- data.frame(x = x, y = sin(2 * x) + s * e0)[1:40, ]
  n <- 40L
  reml <- function(p) {
    V <- exp(2 * p[1]) * exp(-outer(d$x, d$x, "-")^2 / (2 * exp(2 * p[2])))
    diag(V) <- diag(V) * (1 + 1e-12) + exp(2 * p[3])
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e12)
    one <- backsolve(R, rep(1, n), transpose = TRUE)
    ry <- backsolve(R, d$y, transpose = TRUE)
    xvx <- sum(one^2)
    bh <- sum(one * ry) / xvx
    0.5 * sum((ry - bh * one)^2) + sum(log(diag(R))) + 0.5 * log(xvx) + 0.5 * (n - 1) * log(2 * pi)
  }
  w <- character()
  f <- withCallingHandlers(frm(bf(y ~ gp(x)), data = d, family = gaussian(), REML = TRUE),
                           warning = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 100)); invokeRestart("muffleWarning") })
  p0 <- c(f$estimates$theta, log(sigma(f)[1]))
  cf <- nlminb(p0, reml, control = list(iter.max = 800, eval.max = 1000, rel.tol = 1e-12))
  g <- f$obj$gr(f$opt$par)
  cat(sprintf("NCGPR s=%.0e frm %.5f conv %d msg %s th %s sigma %.5g | cf %.5f th %s sigma %.5g | gap %.2e | grad %s\n",
              s, f$opt$objective, f$opt$convergence, f$opt$message,
              paste(sprintf("%.4f", f$estimates$theta), collapse = ","), sigma(f)[1], cf$objective,
              paste(sprintf("%.4f", cf$par[1:2]), collapse = ","), exp(cf$par[3]), f$opt$objective - cf$objective,
              paste(format(g, digits = 2), collapse = " ")))
  print(unique(w))
}

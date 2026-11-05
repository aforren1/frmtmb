# Two fields and REML at low noise against the closed form, the construction that
# warned in test-gp-noncentered.R. Seed 20261102.
#   Rscript dev/ncgp-two.R <lib or base>
a <- commandArgs(TRUE)
libs <- c("C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (a[1] != "base") libs <- c(a[1], libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261102)
n <- 60L
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
d <- data.frame(x = x, y = sin(2 * x) + 1e-2 * e0)[1:40, ]
n <- 40L
d$z <- round(runif(40, 0, 5), 2)
d$y2 <- d$y + cos(d$z)
marg <- function(p, y) {
  K1 <- exp(2 * p[2]) * exp(-outer(d$x, d$x, "-")^2 / (2 * exp(2 * p[3])))
  K2 <- exp(2 * p[4]) * exp(-outer(d$z, d$z, "-")^2 / (2 * exp(2 * p[5])))
  V <- K1 + K2
  diag(V) <- diag(V) + exp(2 * p[6])
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(1e12)
  r <- backsolve(R, y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
}
for (yn in c("y", "y2")) {
  w <- character()
  f <- withCallingHandlers(frm(as.formula(paste(yn, "~ gp(x) + gp(z)")), data = d, family = gaussian()),
                           warning = function(c) { w <<- c(w, substr(conditionMessage(c), 1, 80)); invokeRestart("muffleWarning") })
  y <- d[[yn]]
  p0 <- c(f$estimates$beta[[1]], f$estimates$theta, log(sigma(f)[1]))
  best <- nlminb(p0, marg, y = y, control = list(iter.max = 800, eval.max = 1000, rel.tol = 1e-12))
  for (st in list(c(0, 0, 0, -1, 0, log(0.05)), c(0, 0, 0, 0, 0, log(0.05)), c(0, 0, 0, -3, 0, log(0.05)))) {
    o <- nlminb(st, marg, y = y, control = list(iter.max = 800, eval.max = 1000, rel.tol = 1e-12))
    if (o$objective < best$objective) best <- o
  }
  cat(sprintf("NCGP2 %s %s frm ll %.4f conv %d th %s sigma %.4g | cf ll %.4f th %s sigma %.4g | gap %.2e | frm grad max %.3g\n",
              yn, a[1], -f$opt$objective, f$opt$convergence, paste(sprintf("%.3f", f$estimates$theta), collapse = ","),
              sigma(f)[1], -best$objective, paste(sprintf("%.3f", best$par[2:5]), collapse = ","), exp(best$par[6]),
              -f$opt$objective + best$objective, max(abs(f$obj$gr(f$opt$par)))))
  print(unique(w))
}

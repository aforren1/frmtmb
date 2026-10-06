# What does the separation fit of test-se-check.R:376 do on this BLAS?
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
cat("BLAS probe:", system.time({m <- matrix(1, 600, 600); m %*% m})[[3]],
    "\n")
set.seed(514)
d <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
d$yb <- as.integer(d$z > 0)
w <- list()
ctl <- frmtmb_control()
if (length(args) > 1) ctl <- {
  m <- as.integer(args[2])
  frmtmb_control(optCtrl = list(eval.max = m, iter.max = m))
}
fit <- withCallingHandlers(
  frm(yb ~ z + x, family = bernoulli(), data = d, control = ctl),
  warning = function(x) {
    w[[length(w) + 1L]] <<- conditionMessage(x)
    invokeRestart("muffleWarning")
  })
cat("warnings:\n"); for (x in w) cat(" -", substr(x, 1, 300), "\n")
cat("opt$par:", format(fit$opt$par, digits = 6), "\n")
cat("convergence:", fit$opt$convergence, " evals:", fit$opt$evals, "\n")
cat("max |grad|:", format(max(abs(fit$obj$gr(fit$opt$par))), digits = 4),
    "\n")
cat("se_explained:", format(fit$cache$se_explained), "\n")
sdr <- frmtmb:::sdr_of(fit)
print(sdr$se_lost)
Hx <- fit$obj$he(fit$opt$par)
print(signif(Hx, 6))
print(eigen((Hx + t(Hx)) / 2, only.values = TRUE)$values)

# Reviewer: test-se-check.R's separation fit (seed 514) at the default
# budget and at 4000, printing the code, evaluations and warnings, and
# what glm() says. Usage: Rscript dev/cifixrev-sep.R <lib|base>
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
cat("BLAS probe (s):",
    system.time({m <- matrix(1, 600, 600); m %*% m})[["elapsed"]], "\n")
warns <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}
set.seed(514)
d <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
d$yb <- as.integer(d$z > 0)
for (b in c(NA, 4000)) {
  ctl <- if (is.na(b)) frmtmb_control() else
    frmtmb_control(optCtrl = list(eval.max = b, iter.max = b))
  r <- warns(frm(yb ~ z + x, family = bernoulli(), data = d, control = ctl))
  o <- r$value$opt
  cat("\nbudget", b, "code", o$convergence, "evals", o$evals,
      "iterations", o$iterations, "message", o$message, "\n")
  cat("beta", format(signif(r$value$estimates$beta, 4)), "\n")
  for (w in r$warnings) cat("  W:", substr(w, 1, 300), "\n")
}
g <- warns(glm(yb ~ z + x, family = binomial(), data = d))
cat("\nglm warnings:", g$warnings, " iter:", g$value$iter, "\n")

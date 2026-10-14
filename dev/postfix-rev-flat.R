# Reviewer of lane postfix: a nonlinear fit whose first attempt fails
# and which an EARLIER recovery start (the placed intercepts with the
# slopes at zero) rescues. The lane computes the scan before any
# restart runs, and every scan evaluation of obj$fn moves the tape's
# last.par.best, from which TMB starts the inner problem. Is the
# rescued fit the same as on 0.70.0?
# Usage: Rscript dev/postfix-rev-flat.R <lib or "base">
a <- commandArgs(TRUE)
lib <- if (length(a) && a[1] != "base") a[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(108)
G <- 30
g <- factor(rep(seq_len(G), each = 20))
u <- rnorm(G, 0, 0.5)
x8 <- runif(G * 20, 0, 2)
z <- rnorm(G * 20)
mu8 <- 5 + 0.3 * z + u[g] + 2 * x8
d8 <- data.frame(x = x8, z = z, g = g,
                 y = rgamma(G * 20, shape = 10, rate = 10 / pmax(mu8, 0.5)))
f8 <- bf(y ~ a + b * x, a ~ 1 + z + (1 | g), b ~ 1, nl = TRUE)
for (k in 1:3) {
  msgs <- character()
  fit <- withCallingHandlers(
    frm(f8, data = d8, family = Gamma("identity"),
        prior = prior(normal(2, 2), nlpar = "a"), verbose = (k == 1)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  if (k == 1) for (m in msgs) cat("  [msg]", substr(m, 1, 200), "\n")
  cat(sprintf("rep %d conv %d logLik %.17g\n", k, fit$opt$convergence,
              as.numeric(logLik(fit))))
  cat("  par:", sprintf("%.17g", fit$opt$par), "\n")
  cat("  iterations", fit$opt$iterations, " evaluations",
      paste(fit$opt$evaluations, collapse = " "), "\n")
}

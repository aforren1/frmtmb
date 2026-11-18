# rel084-kt copy of dev/reviews/squant-review-tape4.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Review of lane squant: does the wrong lF derivative of quantile_elf()'s
# cdf_log (squant-review-tape3.txt) reach a fit? Left-censored rows far
# below the bulk ("< 50" non-detects in data near 100), so that at the
# fit their t = (50 - mu) / (lambda sigma) lies in -690 < t, F < 1e-16.
# Run (worktree root):
#   Rscript dev/reviews/squant-review-tape4.R > dev/reviews/squant-review-tape4.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
fd_grad <- function(fn, par) {
  vapply(seq_along(par), function(k) {
    h <- 1e-6 * max(1, abs(par[k]))
    u <- par
    d <- par
    u[k] <- u[k] + h
    d[k] <- d[k] - h
    (fn(u) - fn(d)) / (2 * h)
  }, 0)
}
for (L in c(95, 90, 80, 60)) {
  for (p in c(0.1, 0.5)) {
    for (lam in c(0.25, 1)) {
      set.seed(20261118)
      n <- 200
      d <- data.frame(x = stats::runif(n))
      d$y <- 100 + 2 * d$x + stats::rnorm(n)
      d$c <- 0
      d$c[1:6] <- -1
      d$y[1:6] <- L
      ws <- character()
      fit <- withCallingHandlers(
        frm(bf(y | cens(c) ~ x, quantile = p),
            family = quantile_elf(smoothing = lam), data = d),
        warning = function(w) {
          ws <<- c(ws, substr(conditionMessage(w), 1, 80))
          invokeRestart("muffleWarning")
        })
      par <- fit$opt$par
      dpv <- lapply(ns$eval_dpars(fit)[[1L]], as.numeric)
      tt <- (L - dpv$mu[1:6]) / (lam * dpv$sigma[1:6])
      gt <- as.numeric(fit$obj$gr(par))
      gf <- fd_grad(fit$obj$fn, par)
      o <- stats::optim(par, fit$obj$fn, function(q) fd_grad(fit$obj$fn, q),
                        method = "BFGS",
                        control = list(reltol = 1e-15, maxit = 1000))
      cat(sprintf(paste0("L %g p %.1f lam %.2f: conv %d warn %d | t of the",
                         " censored rows %.0f to %.0f | max |tape grad| %.3g,",
                         " max |fd grad| %.3g | logLik %.4f, BFGS on fn %.4f |",
                         " par %s -> %s %s\n"),
                  L, p, lam, fit$opt$convergence, length(ws), min(tt),
                  max(tt), max(abs(gt)), max(abs(gf)), -fit$opt$objective,
                  -o$value, paste(sprintf("%.4f", par), collapse = " "),
                  paste(sprintf("%.4f", o$par), collapse = " "),
                  paste(unique(ws), collapse = " | ")))
    }
  }
}
cat("DONE\n")

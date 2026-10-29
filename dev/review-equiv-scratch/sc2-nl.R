# Final check: Gamma and inverse.gaussian nl under their DEFAULT links,
# from default starts and from mu-positive starts, on my own data (seeds
# 5101 and 5102, n = 300). Honored: the nl logLik equals the linear one.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ll <- function(expr) tryCatch(as.numeric(logLik(suppressWarnings(
  suppressMessages(expr)))), error = function(e)
    structure(NA_real_, err = substr(conditionMessage(e), 1, 70)))
for (seed in c(5101, 5102)) for (f in c("Gamma", "inverse.gaussian")) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3))
  d$y <- exp(0.4 + 0.3 * d$x + rnorm(n, 0, 0.4))
  g <- get(f, envir = asNamespace("stats"))
  e0 <- if (f == "Gamma") 1 / mean(d$y) else 1 / mean(d$y)^2
  st <- list(beta = c("b0_(Intercept)" = e0, "b1_(Intercept)" = 0))
  fm <- bf(y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)
  a <- ll(frm(y ~ x, family = g(), data = d))
  b0 <- ll(frm(fm, family = g(), data = d))
  b1 <- ll(frm(fm, family = g(), data = d, start = st))
  cat(sprintf("%-17s seed %d  linear %.4f  nl default start %s  nl positive start %.4f  equal %s\n",
              f, seed, a, if (is.na(b0)) paste("FAILS:", attr(b0, "err"))
              else sprintf("%.4f", b0), b1,
              isTRUE(abs(a - b1) < 1e-4 * abs(a))))
}

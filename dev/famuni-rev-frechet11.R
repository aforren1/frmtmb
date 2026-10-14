# Reviewer probe: frechet() at nu = 1.1, the seeds whose SEs are NaN.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
nu0 <- 1.1
for (s in 1:60) {
  set.seed(9000 + s)
  d <- data.frame(x = rnorm(300))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(300)^(-1 / nu0)
  ws <- character(0); ms <- character(0)
  f <- withCallingHandlers(frm(y ~ x, family = frechet(), data = d),
    warning = function(e) { ws <<- c(ws, conditionMessage(e)); invokeRestart("muffleWarning") },
    message = function(e) { ms <<- c(ms, conditionMessage(e)); invokeRestart("muffleMessage") })
  cf <- fixef(f)
  sm <- summary(f)$spec_pars
  if (!all(is.finite(c(cf[, 2], sm[, 2]))) || length(ws)) {
    cat(sprintf("seed %d: conv %d logLik %.4f Intercept %.4f (se %.4g) x %.4f (se %.4g) nu %.5f (se %.4g) link-scale nu par %.4f\n",
                s, f$opt$convergence, as.numeric(logLik(f)), cf[1, 1], cf[1, 2],
                cf[2, 1], cf[2, 2], sm["nu", 1], sm["nu", 2],
                f$opt$par[length(f$opt$par)]))
    cat("   max y", max(d$y), " grad", signif(f$obj$gr(f$opt$par), 3), "\n")
    if (length(ws)) cat("   WARN:", substr(ws, 1, 300), sep = "\n   ")
    if (length(ms)) cat("   MSG:", substr(ms, 1, 400), sep = "\n   ")
    cat("\n")
  }
}

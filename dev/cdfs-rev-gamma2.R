# Reviewer: the NaN points of gamma_cdf_log() at large shapes, and
# censored Gamma fits at large shapes against an independent ML fit.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (a in c(200, 1e3, 1e4, 3e4, 1e5, 1e6)) {
  r <- c(1e-12, 1e-9, 1e-6, 1e-4, 1e-3, 0.01, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9,
         1.5, 2, 5, 30, 1000)
  p <- ns$gamma_cdf_log(r, 1, a)
  rF <- stats::pgamma(r * a, a, log.p = TRUE)
  rS <- stats::pgamma(r * a, a, lower.tail = FALSE, log.p = TRUE)
  bad <- !is.finite(p$lF) | !is.finite(p$lS) |
    abs(p$lF - rF) / pmax(1, abs(rF)) > 1e-10 |
    abs(p$lS - rS) / pmax(1, abs(rS)) > 1e-10
  cat(sprintf("shape %g: bad at q/mu = %s\n", a,
              paste(r[bad], collapse = ", ")))
  if (any(bad)) print(data.frame(q_over_mu = r[bad], lF = p$lF[bad],
                                 ref_lF = rF[bad], lS = p$lS[bad],
                                 ref_lS = rS[bad]), digits = 8)
}
cat("== smallest shape with a NaN at q/mu = 0.5, 0.1, 1e-3\n")
for (r in c(0.5, 0.1, 1e-3)) {
  sh <- 10^seq(1, 7, by = 0.05)
  v <- vapply(sh, function(a) ns$gamma_cdf_log(r, 1, a)$lF, 0)
  rf <- stats::pgamma(r * sh, sh, log.p = TRUE)
  e <- abs(v - rf) / pmax(1, abs(rf))
  i <- which(!is.finite(v) | e > 1e-10)
  cat(sprintf("q/mu %g: first bad shape %s (of %d bad)\n", r,
              if (length(i)) format(sh[i[1]]) else "none", length(i)))
}

cat("== censored Gamma fits vs an independent stats::pgamma ML fit\n")
ind_fit <- function(d) {
  nll <- function(p) {
    mu <- exp(p[1] + p[2] * d$x)
    sh <- exp(p[3])
    ll <- ifelse(d$cc == 0, stats::dgamma(d$y, sh, sh / mu, log = TRUE),
                 stats::pgamma(d$y, sh, sh / mu, lower.tail = FALSE,
                               log.p = TRUE))
    -sum(ll)
  }
  o <- stats::optim(c(log(mean(d$y)), 0, log(10)), nll, method = "BFGS",
                    control = list(maxit = 2000, reltol = 1e-14))
  o2 <- stats::nlminb(o$par, nll)
  c(o2$par, -o2$objective)
}
for (shp in c(5, 300, 3e3, 3e4, 1e5)) {
  set.seed(3)
  d <- data.frame(x = rnorm(200))
  d$y <- rgamma(200, shape = shp, rate = shp / exp(1 + 0.1 * d$x))
  qq <- quantile(d$y, 0.8)
  d$cc <- as.integer(d$y > qq)
  d$y[d$cc == 1] <- qq
  r <- tryCatch({
    fit <- withCallingHandlers(
      frm(y | cens(cc) ~ x, family = Gamma(link = "log"), data = d),
      warning = function(w) {
        cat("  WARN:", substr(conditionMessage(w), 1, 90), "\n")
        invokeRestart("muffleWarning")
      })
    c(fit$opt$par[1:3], as.numeric(logLik(fit)), fit$opt$convergence)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n"); rep(NA, 5)
  })
  ii <- ind_fit(d)
  cat(sprintf("shape %g: frmtmb par %s ll %.8f conv %s | indep par %s ll %.8f\n",
              shp, paste(signif(r[1:3], 8), collapse = " "), r[4], r[5],
              paste(signif(ii[1:3], 8), collapse = " "), ii[4]))
}

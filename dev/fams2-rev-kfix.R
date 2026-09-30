# Reviewer: with kappa held at a constant, the interior rows' dbeta()
# reads a constant first argument, so the fit can reach high phi; there
# the boundary rows go through log_ibeta_half()'s 50-step fraction near
# its switch. Does the fit's logLik and optimum move? Seed 4200 data as
# in dev/fams2-rev-dens.R section B, mu 0.047, kappa 0.05.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ref_xbeta <- function(y, mu, phi, kappa) {
  a <- mu * phi; b <- (1 - mu) * phi; d <- 1 + 2 * kappa
  ifelse(y <= 0, stats::pbeta(kappa / d, a, b, log.p = TRUE),
         ifelse(y >= 1,
                stats::pbeta((1 + kappa) / d, a, b, lower.tail = FALSE,
                             log.p = TRUE),
                stats::dbeta((y + kappa) / d, a, b, log = TRUE) - log(d)))
}
for (cfg in list(c(2e3, 0.047), c(1e4, 0.047), c(2e4, 0.047), c(5e4, 0.047), c(2e5, 0.047), c(1e4, 0.05 / 1.1), c(2e4, 0.05 / 1.1), c(5e4, 0.05 / 1.1), c(2e5, 0.05 / 1.1))) {
  phi <- cfg[1]
  set.seed(4200)
  n <- 1000; mu <- cfg[2]; kap <- 0.05
  z <- rbeta(n, mu * phi, (1 - mu) * phi)
  y <- pmin(pmax((1 + 2 * kap) * z - kap, 0), 1)
  f <- tryCatch(frm(bf(y ~ 1, kappa = 0.05), family = xbeta(),
                    data = data.frame(y = y)),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("phi %.0e: ERROR %s\n", phi, substr(conditionMessage(f), 1, 80)))
    next
  }
  dp <- frmtmb:::eval_dpars(f)[[1]]
  m <- dp$mu[1]; ph <- dp$phi[1]
  ll <- as.numeric(logLik(f))
  ref <- sum(ref_xbeta(y, m, ph, kap))
  nll <- function(p) -sum(ref_xbeta(y, plogis(p[1]), exp(p[2]), kap))
  o <- optim(c(qlogis(m), log(ph)), nll, method = "Nelder-Mead",
             control = list(reltol = 1e-15, maxit = 5000))
  # the boundary-row term alone
  a <- m * ph; b <- (1 - m) * ph; q <- kap / (1 + 2 * kap)
  cat(sprintf(paste0("true mu %.5f ", "phi %.0e: zeros %d; fit mu %.6f phi %.1f | logLik %.8f, ",
                     "reference at the estimates %.8f, diff %.2e | per-row ",
                     "P(0) error %.2e | reference optimum %.8f at mu %.6f ",
                     "phi %.1f\n"),
              mu, phi, sum(y == 0), m, ph, ll, ref, ll - ref,
              frmtmb:::log_ibeta_half(q, a, b) - pbeta(q, a, b, log.p = TRUE),
              -o$value, plogis(o$par[1]), exp(o$par[2])))
}

# Reviewer: log_igamma_two() is NaN below the mode once a > ~12600
# (0 * Inf in a blend arm weighted zero). Who reaches it: Gamma and
# hurdle_gamma (shape a), and the poisson / zero-inflated / hurdle
# poisson pair (a = q + 1, x = mu). Run with ARM=lane or ARM=base.
arm <- Sys.getenv("ARM", "lane")
lp <- c("C:/Users/adf44/source/r/rellib-r10",
        "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") lp <- c("C:/Users/adf44/source/r/wt-cdfs-rev-lib", lp)
.libPaths(lp)
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
if (arm == "lane") {
  cat("== NaN map of log_igamma_two(x, a): count of NaN lP or lQ\n")
  for (a in c(5e3, 1e4, 1.2e4, 1.3e4, 2e4, 5e4, 1e5, 1e6)) {
    r <- c(0.01, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 1, 1.02,
           1.05, 1.2, 1.5, 3)
    p <- ns$log_igamma_two(r * a, a)
    bad <- r[is.nan(p$lP) | is.nan(p$lQ)]
    cat(sprintf("a %-7g NaN at x/a = %s\n", a, paste(bad, collapse = " ")))
  }
  cat("== the gl arm alone at a = 1e5, x = 0.7 (a + 1): lQ\n")
  print(ns$igamma_gl_tails(0.7 * (1e5 + 1), 1e5))
}
cat("== poisson lccdf slot (P(Y > q)), vs stats::ppois\n")
fp <- ns$as_frmtmb_family(poisson())
for (pr in list(c(1e5, 7e4), c(1e5, 9e4), c(3e4, 2e4), c(2e4, 1.5e4),
                c(1e5, 1e5 - 600))) {
  q <- pr[1]; mu <- pr[2]
  ref <- stats::ppois(q, mu, lower.tail = FALSE, log.p = TRUE)
  got <- if (!is.null(fp$lccdf)) fp$lccdf(q, list(mu = mu), list()) else
    log1p(-fp$lcdf(q, list(mu = mu), list()))
  cat(sprintf("q %g mu %g: slot %s  ref %.10g\n", q, mu, format(got, digits = 10),
              ref))
}
cat("== a right-censored poisson fit with large counts\n")
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n))
d$y <- rpois(n, exp(11 + 0.3 * d$x))
d$cc <- 0
# a detector that saturates: rows above 75000 are recorded as 75000
d$cc[d$y > 75000] <- 1
d$y[d$cc == 1] <- 75000
cat("rows censored:", sum(d$cc), "of", n, "\n")
r <- tryCatch({
  fit <- withCallingHandlers(
    frm(y | cens(cc) ~ x, family = poisson(), data = d),
    warning = function(w) {
      cat("  WARN:", substr(conditionMessage(w), 1, 100), "\n")
      invokeRestart("muffleWarning")
    })
  cat(sprintf("coef %s  logLik %.6f conv %d\n",
              paste(signif(fit$opt$par, 10), collapse = " "),
              as.numeric(logLik(fit)), fit$opt$convergence))
  cat("objective at the estimate:", fit$obj$fn(fit$opt$par), " gradient:",
      fit$obj$gr(fit$opt$par), "\n")
}, error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
# independent ML: P(Y >= 75000) with the inclusive convention frmtmb uses
nll <- function(p) {
  mu <- exp(p[1] + p[2] * d$x)
  -sum(ifelse(d$cc == 0, dpois(d$y, mu, log = TRUE),
              ppois(d$y - 1, mu, lower.tail = FALSE, log.p = TRUE)))
}
o <- nlminb(c(11, 0.3), nll)
cat(sprintf("independent: coef %s logLik %.6f\n",
            paste(signif(o$par, 10), collapse = " "), -o$objective))

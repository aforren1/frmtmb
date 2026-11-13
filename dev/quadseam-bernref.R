# dev/quadseam-acc.R's bernoulli design, seeds 1..10: each build's
# log-likelihood at its own estimate against an independent reference
# there (stats::integrate() per group, plain R, rel.tol 1e-12), so the
# rule's own error is read apart from where the optimizer stopped.
#   Rscript dev/quadseam-bernref.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
for (seed in 1:10) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:100, each = 4)), x = rnorm(400))
  d$y <- rbinom(400, 1, plogis(-0.5 + 0.7 * d$x + rnorm(100, 0, 1.2)[d$g]))
  f <- suppressMessages(frm(y ~ x + (1 | g), family = bernoulli(),
                            data = d, quadrature = TRUE))
  b <- f$estimates$beta
  s <- exp(f$estimates$theta[1])
  ref <- sum(vapply(split(d, d$g), function(dg) {
    h <- function(u) vapply(u, function(ui) {
      exp(sum(dbinom(dg$y, 1, plogis(b[1] + b[2] * dg$x + ui), log = TRUE)) +
            dnorm(ui, 0, s, log = TRUE))
    }, 0)
    log(integrate(h, -Inf, Inf, rel.tol = 1e-12)$value)
  }, 0))
  cat(sprintf("SEED %2d fit-ref %10.3e\n", seed,
              as.numeric(logLik(f)) - ref))
}

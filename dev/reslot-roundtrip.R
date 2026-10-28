# Lane reslot, item 2: the sampling table of a residual-correlation term
# prints brms's (flat) with lb/ub, and the table passed back as prior =
# gives the same sampled density as the defaults. Compared through the
# retaped objective and its gradient at 300 points per model (seed 7),
# never through draws (dev/rtmb-pitfalls.md item 23).
# Usage: Rscript dev/reslot-roundtrip.R before|after
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages({
  library(frmtmb.sample)
  library(frmtmb.spatial)
})
set.seed(4)
ng <- 12L
k <- 5L
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                x = rnorm(ng * k))
d$y <- 1 + 0.5 * d$x + rnorm(ng)[d$g] * 0.3 + rnorm(ng * k)
W <- matrix(0, 10, 10)
for (i in 1:9) W[i, i + 1] <- W[i + 1, i] <- 1
ds <- data.frame(x = rnorm(10))
ds$y <- 1 + ds$x + rnorm(10)
cases <- list(
  ar_cov = list(f = y ~ x + ar(t, g, cov = TRUE), d = d),
  ma_cov = list(f = y ~ x + ma(t, g, cov = TRUE), d = d),
  arma_cov = list(f = y ~ x + arma(t, g, cov = TRUE), d = d),
  ar2_cov = list(f = y ~ x + ar(t, g, p = 2, cov = TRUE), d = d),
  cosy = list(f = y ~ x + cosy(t, g), d = d),
  sar = list(f = y ~ x + sar(W), d = ds, d2 = list(W = W)))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  d2 <- cs$d2 %||% list()
  cat("\n==", nm, "\n")
  tab <- default_prior(cs$f, data = cs$d, data2 = d2, route = "sample")
  ac <- tab[!tab$class %in% c("b", "Intercept", "sigma"), ]
  print(as.data.frame(ac)[, c("prior", "class", "lb", "ub")],
        row.names = FALSE)
  fit <- suppressWarnings(suppressMessages(
    frm(cs$f, data = cs$d, data2 = d2)))
  r0 <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))
  cat("prior_summary of the defaults:\n")
  pe <- as.data.frame(r0$effective)
  print(pe[!pe$class %in% c("Intercept", "sigma"),
           c("prior", "class", "lb", "ub")], row.names = FALSE)
  r1 <- tryCatch(suppressMessages(
    frmtmb.sample:::sample_resolve_priors(fit, tab)),
    error = function(e) e)
  if (inherits(r1, "error")) {
    cat("table as prior = ERROR:", conditionMessage(r1), "\n")
    next
  }
  ch0 <- autocor_bound_chart(fit, r0$ri$lower, r0$ri$upper)
  ch1 <- autocor_bound_chart(fit, r1$ri$lower, r1$ri$upper)
  b0 <- frmtmb:::resolve_bounds(fit, r0$ri$lower, r0$ri$upper)
  b1 <- frmtmb:::resolve_bounds(fit, r1$ri$lower, r1$ri$upper)
  cat("bounds identical:", identical(b0, b1), " charts:", length(ch0),
      length(ch1), "\n")
  o0 <- frmtmb.sample:::prior_augmented_obj(fit, r0$ri$entries, ch0)
  o1 <- frmtmb.sample:::prior_augmented_obj(fit, r1$ri$entries, ch1)
  set.seed(7)
  p0 <- o0$par
  dfn <- 0
  dgr <- 0
  nfin <- 0
  for (j in seq_len(300)) {
    p <- p0 + rnorm(length(p0), 0, 0.5)
    f0 <- o0$fn(p)
    f1 <- o1$fn(p)
    if (is.finite(f0)) nfin <- nfin + 1
    dfn <- max(dfn, abs(f0 - f1) / max(1, abs(f0)), na.rm = TRUE)
    dgr <- max(dgr, max(abs(o0$gr(p) - o1$gr(p))) /
                 max(1, max(abs(o0$gr(p)))), na.rm = TRUE)
  }
  cat(sprintf("300 points (%d finite): max rel |fn diff| %.3g, max rel |gr diff| %.3g\n",
              nfin, dfn, dgr))
}

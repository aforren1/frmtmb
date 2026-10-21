# Release 0.73.0 consolidation: copy of lane adhess's battery script, with the "lane" arm on wt-release-lib (the merged build) and "base" on rellib-r10 (0.72.0); paths point at the release tree.
# Lane adhess, copied from lane illse2's copy (from dev/illse-sf-rev4-bitwise.R): copy of dev/rel070-seflags-rev4-bitwise.R by dev/illse-sf-make.R (lane arm: the lane library over
# rellib-r10; base arm: rellib-r10).
# Reviewer of lane seflags, final check of the second defect's fix: on
# fits where the check re-measured no span (nothing lost, tier 1
# stands), predictions must still read sdreport()'s inverse, bit for
# bit. Per fit: whether the lane kept a repaired covariance
# (sdr$cov_fixed_prop), and fixef() SEs, vcov(), fitted(newdata) and
# predict(newdata) SEs, saved for comparison with the other arm.
#   Rscript dev/seflags-rev4-bitwise.R <lane|base|op3> <out.rds>
source("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel073-sf-rev4-common.R")
out <- commandArgs(trailingOnly = TRUE)[2]
mk <- function(s, n = 300) {
  set.seed(s)
  d <- data.frame(x = stats::rnorm(n), z = stats::runif(n),
                  g = factor(rep(1:20, length.out = n)))
  u <- stats::rnorm(20, 0, 0.5)[d$g]
  us <- stats::rnorm(20, 0, 0.3)[d$g]
  d$y <- 1 + 0.5 * d$x + sin(3 * d$z) + u + us * d$x + stats::rnorm(n, 0, 0.4)
  d$yb <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * d$x + u))
  d$yc <- stats::rpois(n, exp(0.2 + 0.4 * d$x + u))
  d$yg <- stats::rgamma(n, 4, 4 / exp(0.5 + 0.3 * d$x))
  d$yn <- stats::rnbinom(n, mu = exp(1 + 0.3 * d$x + u), size = 2)
  d
}
specs <- list(
  lm = list(y ~ x + z, gaussian()),
  ri = list(y ~ x + z + (1 | g), gaussian()),
  rs = list(y ~ x + z + (1 + x | g), gaussian()),
  smooth = list(y ~ x + s(z) + (1 | g), gaussian()),
  sigma = list(bf(y ~ x + (1 | g), sigma ~ z), gaussian()),
  bern = list(yb ~ x + z + (1 | g), bernoulli()),
  pois = list(yc ~ x + z + (1 | g), poisson()),
  gamma = list(yg ~ x + z, Gamma(link = "log")),
  negbin = list(yn ~ x + (1 | g), negbinomial()),
  poly3 = list(y ~ poly(x, 3) + z + (1 | g), gaussian()))
res <- list()
for (s in 1:3) {
  d <- mk(s)
  nd <- d[c(1, 7, 50, 120), ]
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    f <- suppressWarnings(suppressMessages(frm(sp[[1]], family = sp[[2]],
                                               data = d)))
    sdr <- suppressWarnings(ns$sdr_of(f))
    g <- function(e) tryCatch(suppressWarnings(e), error = function(e) NULL)
    res[[paste(nm, s)]] <- list(
      code = f$opt$convergence,
      lost = sdr$se_lost,
      repaired = !is.null(sdr$cov_fixed_prop),
      fixef = g(fixef(f)[, "Est.Error"]),
      vcov = g(stats::vcov(f)),
      fitted = g(fitted(f, newdata = nd)[, "Est.Error"]),
      fitted_pop = g(fitted(f, newdata = nd, re_formula = NA)[, "Est.Error"]),
      predict_lin = g(fitted(f, newdata = nd, scale = "linear")[, "Est.Error"]))
    cat(sprintf("%-8s seed %d code %d lost %-3d repaired %s\n", nm, s,
                f$opt$convergence, length(sdr$se_lost),
                !is.null(sdr$cov_fixed_prop)))
  }
}
saveRDS(res, out)

# Lane osa: residuals(type = "osa") on the discrete core families, against
# the exact randomized quantile residual of the fitted CDF.
#
#   Rscript dev/osa-probe.R before   # rellib-r18 only (frmtmb 0.79.1)
#   Rscript dev/osa-probe.R after    # wt-osa-lib first
#
# With no random effect the one-step predictive of row k is the fitted
# row distribution itself, so the OSA residual must lie in
# [qnorm(F(y - 1)), qnorm(F(y))], and with TMB's own uniform draws
# (set.seed(123); runif(n), oneStepPredict()'s default seed) it must equal
# qnorm(F(y - 1) + (1 - U) * p(y)) up to rounding. The second is printed
# as the largest absolute difference; the first is the count of rows
# outside the interval. Data seed 9, n = 150, 10 trials (the 0.79.1
# record's construction for binomial); the other designs reuse it.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "after"
lib <- c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
         "C:/Users/adf44/source/r/rellib-r18",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressMessages(library(frmtmb))
cat("arm:", arm, " frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), " RTMBdist",
    format(packageVersion("RTMBdist")), "\n")

pbb <- function(q, size, a, b) {
  vapply(seq_along(q), function(i) {
    if (q[i] < 0) return(0)
    k <- 0:q[i]
    sum(exp(lchoose(size[i], k) + lbeta(k + a[i], size[i] - k + b[i]) -
              lbeta(a[i], b[i])))
  }, 0)
}

truth_cdf <- function(fam, y, dp, size) {
  switch(fam,
    poisson = list(lo = ppois(y - 1, dp$mu), hi = ppois(y, dp$mu)),
    binomial = list(lo = pbinom(y - 1, size, dp$mu),
                    hi = pbinom(y, size, dp$mu)),
    bernoulli = list(lo = pbinom(y - 1, 1, dp$mu), hi = pbinom(y, 1, dp$mu)),
    negbinomial = list(lo = pnbinom(y - 1, size = dp$shape, mu = dp$mu),
                       hi = pnbinom(y, size = dp$shape, mu = dp$mu)),
    geometric = list(lo = pnbinom(y - 1, size = 1, mu = dp$mu),
                     hi = pnbinom(y, size = 1, mu = dp$mu)),
    beta_binomial = {
      a <- dp$mu * dp$phi
      b <- (1 - dp$mu) * dp$phi
      list(lo = pbb(y - 1, size, a, b), hi = pbb(y, size, a, b))
    })
}

set.seed(9)
n <- 150
x <- rnorm(n)
size <- rep(10L, n)
eta <- -0.3 + 0.8 * x
pb <- plogis(eta)
dat <- list(
  poisson = rpois(n, exp(0.8 + 0.5 * x)),
  binomial = rbinom(n, size, pb),
  bernoulli = rbinom(n, 1, pb),
  negbinomial = rnbinom(n, size = 2, mu = exp(0.8 + 0.5 * x)),
  geometric = rnbinom(n, size = 1, mu = exp(0.8 + 0.5 * x)),
  beta_binomial = rbinom(n, size, rbeta(n, pb * 8, (1 - pb) * 8)))

U <- local({ set.seed(123); runif(n) })
rows <- list()
for (fam in names(dat)) {
  d <- data.frame(y = dat[[fam]], x = x, k = size)
  form <- if (fam %in% c("binomial", "beta_binomial")) {
    bf(y | trials(k) ~ x)
  } else bf(y ~ x)
  famobj <- frmtmb:::family_registry[[fam]]()
  fit <- frm(form, family = famobj, data = d)
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(residuals(fit, type = "osa")[, "Estimate"],
                error = function(e) structure(conditionMessage(e),
                                              class = "err"))
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "err")) {
    cat(sprintf("%-14s ERROR %s\n", fam, substr(r, 1, 100)))
    rows[[fam]] <- data.frame(family = fam, nan = NA, outside = NA,
                              maxdiff = NA, mean = NA, sd = NA, secs = el,
                              note = substr(r, 1, 60))
    next
  }
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  dp <- lapply(dp, as.numeric)
  tc <- truth_cdf(fam, d$y, dp, size)
  qlo <- qnorm(tc$lo)
  qhi <- qnorm(tc$hi)
  fin <- is.finite(r)
  wid <- pmax(qhi - qlo, 0)
  # relative slack: 1e-6 of the row's own interval width
  outside <- sum(fin & (r < qlo - 1e-6 * wid | r > qhi + 1e-6 * wid))
  exact <- qnorm(tc$hi - U * (tc$hi - tc$lo))
  rows[[fam]] <- data.frame(
    family = fam, nan = sum(!fin), outside = outside,
    maxdiff = if (any(fin)) max(abs(r - exact)[fin]) else NA,
    mean = mean(r[fin]), sd = sd(r[fin]), secs = el, note = "")
}
res <- do.call(rbind, rows)
rownames(res) <- NULL
cat("BEGIN osa-probe", arm, "\n")
print(res, digits = 4)
cat("END osa-probe\n")

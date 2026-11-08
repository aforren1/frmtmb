# Lane osa: the edges of residuals(type = "osa") that the audit's fixed
# design does not reach: trunc() on counts, case weights, the explicit
# methods, and the refusals.
#
#   Rscript dev/osa-edge.R before|after
#
# Data seed 9 throughout (the audit's x), n = 150. Truth is the fitted
# row distribution's CDF from stats:: functions, inside the trunc()
# window where there is one; "maxdiff" is against the randomized value
# with TMB's own uniforms (set.seed(123); runif(n)), and "outside" counts
# rows outside [qnorm(F(y - 1)), qnorm(F(y))].
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "after"
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " frmtmb from", dirname(find.package("frmtmb")),
    " RTMBdist", format(packageVersion("RTMBdist")), "\n")
set.seed(9)
n <- 150
x <- rnorm(n)
U <- local({ set.seed(123); runif(n) })
run <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
                           warning = function(c) {
                             w <<- c(w, conditionMessage(c))
                             invokeRestart("muffleWarning")
                           })
  list(v = v, w = w)
}
show <- function(label, res, lo = NULL, hi = NULL) {
  if (inherits(res$v, "error")) {
    cat(sprintf("%-34s ERROR %s\n", label,
                substr(gsub("\n", " ", conditionMessage(res$v)), 1, 150)))
    return(invisible())
  }
  r <- if (is.matrix(res$v)) res$v[, "Estimate"] else res$v
  msg <- sprintf("nan %d", sum(!is.finite(r)))
  if (!is.null(lo)) {
    qlo <- qnorm(lo)
    qhi <- qnorm(pmin(hi, 1))
    wid <- pmax(qhi - qlo, 0)
    wid[!is.finite(wid)] <- 1
    fin <- is.finite(r)
    ex <- qnorm(pmin(hi, 1) - U * (hi - lo))
    ok <- fin & is.finite(ex)
    msg <- paste(msg, sprintf("outside %d maxdiff %.3g", sum(fin & (
      r < qlo - 1e-6 * wid | r > qhi + 1e-6 * wid)),
      max(abs(r - ex)[ok])))
  }
  if (length(res$w)) msg <- paste(msg, "WARN:", substr(res$w[1], 1, 120))
  cat(sprintf("%-34s %s\n", label, msg))
}
osa <- function(f, ...) residuals(f, type = "osa", ...)
dpar <- function(f) lapply(frmtmb:::eval_dpars(f)[[1]], as.numeric)

cat("BEGIN osa-edge", arm, "\n")
# zero-truncated poisson, trunc(lb = 1)
d <- data.frame(x = x, y = qpois(runif(n, dpois(0, exp(1 + 0.4 * x)), 1),
                                 exp(1 + 0.4 * x)))
f <- frm(y | trunc(lb = 1) ~ x, family = poisson(), data = d)
m <- dpar(f)$mu
Z <- 1 - ppois(0, m)
show("poisson trunc(lb = 1)", run(osa(f)),
     (ppois(d$y - 1, m) - ppois(0, m)) / Z, (ppois(d$y, m) - ppois(0, m)) / Z)
# binomial truncated above, trunc(ub = 7) of 10 trials
pb <- plogis(-0.3 + 0.8 * x)
yb <- rbinom(n, 10, pb)
while (any(yb > 7)) yb[yb > 7] <- rbinom(sum(yb > 7), 10, pb[yb > 7])
d <- data.frame(x = x, y = yb, k = 10L)
f <- frm(y | trials(k) + trunc(ub = 7) ~ x, family = binomial(), data = d)
p <- dpar(f)$mu
Z <- pbinom(7, 10, p)
show("binomial trunc(ub = 7)", run(osa(f)),
     pbinom(d$y - 1, 10, p) / Z, pbinom(d$y, 10, p) / Z)
# row-varying lower bound on a count
lb <- rep(c(1, 2), length.out = n)
mu0 <- exp(1.2 + 0.4 * x)
yv <- qpois(runif(n, ppois(lb - 1, mu0), 1), mu0)
d <- data.frame(x = x, y = yv, lb = lb)
f <- frm(y | trunc(lb = lb) ~ x, family = poisson(), data = d)
m <- dpar(f)$mu
Z <- 1 - ppois(lb - 1, m)
show("poisson trunc(lb = row-varying)", run(osa(f)),
     (ppois(d$y - 1, m) - ppois(lb - 1, m)) / Z,
     (ppois(d$y, m) - ppois(lb - 1, m)) / Z)

# case weights, no random effect: the row's own CDF
d <- data.frame(x = x, y = rpois(n, exp(0.8 + 0.5 * x)),
                w = rep(c(1, 3), length.out = n))
f <- frm(y | weights(w) ~ x, family = poisson(), data = d)
m <- dpar(f)$mu
show("poisson weights(1, 3)", run(osa(f)), ppois(d$y - 1, m), ppois(d$y, m))
dg <- data.frame(x = x, y = 1 + 0.5 * x + rnorm(n),
                 w = rep(c(1, 3), length.out = n))
fg <- frm(y | weights(w) ~ x, family = gaussian(), data = dg)
rg <- run(osa(fg))
dpg <- dpar(fg)
z <- (dg$y - dpg$mu) / dpg$sigma
if (!inherits(rg$v, "error")) {
  r <- rg$v[, "Estimate"]
  cat(sprintf("%-34s max|r - z| %.3g  max|r - z sqrt(w)| %.3g\n",
              "gaussian weights(1, 3)", max(abs(r - z)),
              max(abs(r - z * sqrt(dg$w)))))
  fc <- run(osa(fg, osa_method = "cdf"))
  if (!inherits(fc$v, "error")) {
    cat(sprintf("%-34s max|r - z| %.3g\n", "gaussian weights, cdf",
                max(abs(fc$v[, "Estimate"] - z))))
  } else show("gaussian weights, cdf", fc)
}

# explicit methods on counts
d <- data.frame(x = x, y = rbinom(n, 10, pb), k = 10L)
f <- frm(y | trials(k) ~ x, family = binomial(), data = d)
p <- dpar(f)$mu
show("binomial oneStepGeneric 0:10",
     run(osa(f, osa_method = "oneStepGeneric", discreteSupport = 0:10)),
     pbinom(d$y - 1, 10, p), pbinom(d$y, 10, p))
show("binomial oneStepGeneric no support",
     run(osa(f, osa_method = "oneStepGeneric")))
show("binomial fullGaussian", run(osa(f, osa_method = "fullGaussian")))
show("binomial oneStepGaussianOffMode",
     run(osa(f, osa_method = "oneStepGaussianOffMode")))
d$y <- rbinom(n, 10, rbeta(n, pb * 8, (1 - pb) * 8))
f <- frm(y | trials(k) ~ x, family = beta_binomial(), data = d)
dp <- dpar(f)
pbb <- function(q, a, b) vapply(seq_along(q), function(i) {
  if (q[i] < 0) return(0)
  k <- 0:q[i]
  sum(exp(lchoose(10, k) + lbeta(k + a[i], 10 - k + b[i]) - lbeta(a[i], b[i])))
}, 0)
a <- dp$mu * dp$phi
b <- (1 - dp$mu) * dp$phi
show("beta_binomial default", run(osa(f)), pbb(d$y - 1, a, b), pbb(d$y, a, b))
show("beta_binomial oneStepGeneric 0:10",
     run(osa(f, osa_method = "oneStepGeneric", discreteSupport = 0:10)),
     pbb(d$y - 1, a, b), pbb(d$y, a, b))
show("beta_binomial oneStepGaussianOffMode",
     run(osa(f, osa_method = "oneStepGaussianOffMode")))

# a family with no CDF, and the warning on a non-finite residual
dh <- data.frame(x = x, y = 1 + 0.7 * x + rnorm(n))
fh <- frm(y ~ x, family = huber(), data = dh)
show("huber default", run(osa(fh)))
dig <- data.frame(x = x, y = rgamma(n, 3, 3 / exp(0.3 + 0.4 * x)))
fi <- frm(y ~ x, family = inverse.gaussian(), data = dig)
show("inverse.gaussian oneStepGeneric", run(osa(fi, osa_method =
                                                  "oneStepGeneric")))
show("inverse.gaussian default", run(osa(fi)))

# RTMBdist-5: dbetabinom() on the tape as the shapes grow
for (lphi in c(10, 39, 45)) {
  tp <- RTMB::MakeTape(function(p) {
    sum(RTMBdist::dbetabinom(c(5, 3, 7), 10, 0.5 * exp(p), 0.5 * exp(p),
                             log = TRUE))
  }, lphi)
  cat(sprintf("RTMBdist-5 dbetabinom tape  log phi %2d: %.6f   binomial %.6f\n",
              lphi, tp(lphi), sum(dbinom(c(5, 3, 7), 10, 0.5, log = TRUE))))
}
cat("END osa-edge\n")

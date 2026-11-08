# Lane osa: residuals(type = "osa") with a random intercept, against the
# exact one-step predictive.
#
#   Rscript dev/osa-re.R before|after [reps]
#
# Design: 20 groups of 10 rows, group-major, x ~ N(0, 1). A base data
# set (seed 11) is fitted; `reps` responses are then simulated from that
# fit with the group effects redrawn (simulate(re_formula = NA, seed =
# 21)), and each is refitted. For each refit:
#
# - exact: with the parameters held at the refit's estimates, the
#   one-step predictive of row k is int F(y_k | u) p(u | earlier rows of
#   its group) du (groups are independent, so TMB's conditioning on every
#   earlier row reduces to the group's). It is formed here on a grid of
#   801 points over +-8 random-effect sd. With TMB's own uniforms
#   (set.seed(123); runif(n)) the exact residual is
#   qnorm(F(y) - U p(y)); "maxdiff" and "meddiff" compare it with the OSA
#   residual, which is the Laplace approximation of the same integral.
# - calibration: the OSA residuals of all refits pooled, KS against
#   N(0, 1), with mean and sd.
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "after"
reps <- as.integer(commandArgs(TRUE)[2])
if (is.na(reps)) reps <- 10L
.libPaths(c(if (arm == "after") "C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " frmtmb from", dirname(find.package("frmtmb")), "\n")

G <- 20L
K <- 10L
n <- G * K
size <- 10L
U <- local({ set.seed(123); runif(n) })

pbb <- function(q, sz, a, b) {
  vapply(seq_along(q), function(i) {
    if (q[i] < 0) return(0)
    k <- 0:q[i]
    sum(exp(lchoose(sz, k) + lbeta(k + a[i], sz - k + b[i]) - lbeta(a[i], b[i])))
  }, 0)
}
dbb <- function(y, sz, a, b) {
  exp(lchoose(sz, y) + lbeta(y + a, sz - y + b) - lbeta(a, b))
}
# F(y | eta) and p(y | eta) on a vector of linear predictors
fam_fun <- list(
  poisson = list(link = log, F = function(y, eta, dp) ppois(y, exp(eta)),
                 p = function(y, eta, dp) dpois(y, exp(eta))),
  bernoulli = list(link = qlogis,
                   F = function(y, eta, dp) pbinom(y, 1, plogis(eta)),
                   p = function(y, eta, dp) dbinom(y, 1, plogis(eta))),
  binomial = list(link = qlogis,
                  F = function(y, eta, dp) pbinom(y, size, plogis(eta)),
                  p = function(y, eta, dp) dbinom(y, size, plogis(eta))),
  negbinomial = list(link = log,
                     F = function(y, eta, dp) {
                       pnbinom(y, size = dp$shape, mu = exp(eta))
                     },
                     p = function(y, eta, dp) {
                       dnbinom(y, size = dp$shape, mu = exp(eta))
                     }),
  beta_binomial = list(link = qlogis,
                       F = function(y, eta, dp) {
                         m <- plogis(eta)
                         pbb(rep(y, length(eta)), size, m * dp$phi,
                             (1 - m) * dp$phi)
                       },
                       p = function(y, eta, dp) {
                         m <- plogis(eta)
                         dbb(y, size, m * dp$phi, (1 - m) * dp$phi)
                       }))

exact_osa <- function(fit, d, fam) {
  ff <- fam_fun[[fam]]
  dp <- lapply(frmtmb:::eval_dpars(fit)[[1]], as.numeric)
  dps <- list(shape = dp$shape[1], phi = dp$phi[1])
  b <- fit$estimates$b
  eta_fix <- ff$link(dp$mu) - b[as.integer(d$g)]
  sdu <- exp(fit$estimates$theta[[1]])
  u <- seq(-8, 8, length.out = 801) * sdu
  w0 <- dnorm(u, 0, sdu)
  lo <- hi <- numeric(n)
  for (gi in split(seq_len(n), d$g)) {
    w <- w0
    for (i in gi) {
      eta <- eta_fix[i] + u
      p <- ff$p(d$y[i], eta, dps)
      Fh <- ff$F(d$y[i], eta, dps)
      s <- sum(w)
      hi[i] <- sum(w * Fh) / s
      lo[i] <- sum(w * (Fh - p)) / s
      w <- w * p
    }
  }
  qnorm(hi - U * (hi - lo))
}

make_base <- function(fam) {
  set.seed(11)
  d <- data.frame(g = factor(rep(seq_len(G), each = K)), x = rnorm(n),
                  k = size)
  e <- rnorm(G, 0, 0.7)[d$g]
  eta <- 0.2 + 0.5 * d$x + e
  d$y <- switch(fam,
    poisson = rpois(n, exp(eta)),
    bernoulli = rbinom(n, 1, plogis(eta)),
    binomial = rbinom(n, size, plogis(eta)),
    negbinomial = rnbinom(n, size = 3, mu = exp(eta)),
    beta_binomial = {
      m <- plogis(eta)
      rbinom(n, size, rbeta(n, m * 8, (1 - m) * 8))
    })
  d
}
form_of <- function(fam) {
  if (fam %in% c("binomial", "beta_binomial")) {
    bf(y | trials(k) ~ x + (1 | g))
  } else bf(y ~ x + (1 | g))
}

cat("BEGIN osa-re", arm, "reps", reps, "\n")
for (fam in names(fam_fun)) {
  famobj <- frmtmb:::family_registry[[fam]]()
  d0 <- make_base(fam)
  f0 <- suppressWarnings(frm(form_of(fam), family = famobj, data = d0))
  sims <- simulate(f0, nsim = reps, re_formula = NA, seed = 21)
  rr <- ee <- numeric(0)
  nerr <- 0L
  secs <- 0
  for (j in seq_len(reps)) {
    d <- d0
    d$y <- sims[[j]]
    f <- suppressWarnings(frm(form_of(fam), family = famobj, data = d))
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(suppressWarnings(residuals(f, type = "osa")[, "Estimate"]),
                  error = function(e) e)
    secs <- secs + proc.time()[["elapsed"]] - t0
    if (inherits(r, "error")) {
      nerr <- nerr + 1L
      if (nerr == 1L) cat("  ", fam, "ERR", conditionMessage(r), "\n")
      next
    }
    rr <- c(rr, r)
    ee <- c(ee, exact_osa(f, d, fam))
  }
  if (!length(rr)) {
    cat(sprintf("%-14s refused on %d of %d refits\n", fam, nerr, reps))
    next
  }
  ok <- is.finite(rr) & is.finite(ee)
  cat(sprintf(paste("%-14s rows %d nan %d err %d | osa mean %6.3f sd %5.3f",
                    "ks.p %.3g | exact ks.p %.3g | maxdiff %.3g meddiff %.3g",
                    "| %.1fs\n"),
              fam, length(rr), sum(!is.finite(rr)), nerr,
              mean(rr[is.finite(rr)]), sd(rr[is.finite(rr)]),
              suppressWarnings(ks.test(rr[is.finite(rr)], "pnorm")$p.value),
              suppressWarnings(ks.test(ee[is.finite(ee)], "pnorm")$p.value),
              max(abs(rr - ee)[ok]), median(abs(rr - ee)[ok]), secs))
}
cat("END osa-re\n")

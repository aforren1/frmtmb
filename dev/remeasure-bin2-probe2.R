# Lane remeasure, bin 2P groups G3 and G4 (tests.log_lik.R and
# tests.posterior_predict.R through dev/brmsport-bin2-adapters.R): the
# constructions behind the verdicts and the checks that the passes are
# genuine. Seeds are the blocks' own (withr::local_seed(<line>), as the
# generated file sets them).
#
#   Rscript dev/remeasure-bin2-probe2.R \
#     > dev/remeasure-out/probe/bin2-probe2.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(testthat)})
source("dev/brmsport-bin2-adapters.R")
cat("frmtmb", format(packageVersion("frmtmb")), "brms",
    format(packageVersion("brms")), "\n")
inv_logit <- function(x) 1 / (1 + exp(-x))
num <- function(x) formatC(x, digits = 17, format = "g")
cmp <- function(lab, a, b) {
  fin <- is.finite(a) & is.finite(b)
  cat(sprintf(paste("%-26s n=%d finite=%d identical=%s max|diff|=%s",
                    "max rel=%s range(ref)=[%s, %s]\n"),
              lab, length(a), sum(fin), identical(a, b),
              num(max(abs(a - b)[fin], 0)),
              num(max((abs(a - b) / abs(b))[fin], 0)),
              num(min(b)), num(max(b))))
}

cat("\n## 1. The density rows that hold: frmtmb against the reference\n")
cat("(block 'location shift', seed 3; block 'skewed normal', seed 35)\n")
withr::with_seed(3L, {
  ns <- 25
  prep <- structure(list(), class = "brmsprep")
  prep$dpars <- list(mu = matrix(rnorm(ns * 2), ncol = 2),
                     sigma = rchisq(ns, 3), nu = rgamma(ns, 4))
  prep$family <- gaussian()
  prep$data <- list(Y = rnorm(ns))
  ref <- dnorm(prep$data$Y[1], prep$dpars$mu[, 1], prep$dpars$sigma,
               log = TRUE)
  cmp("gaussian (log_lik:19)", brms_adapt_log_lik("gaussian", 1, prep), ref)
  ref_t <- brms::dstudent_t(prep$data$Y[2], prep$dpars$nu,
                            prep$dpars$mu[, 2], prep$dpars$sigma, log = TRUE)
  cmp("student (log_lik:27)", brms_adapt_log_lik("student", 2, prep), ref_t)
  prep$data$weights <- sample(1:10, ns, replace = TRUE)
  cat("weights[1] =", prep$data$weights[1], "\n")
  cmp("weighted (log_lik:32)", brms_adapt_log_lik("gaussian", 1, prep),
      ref * prep$data$weights[1])
})
withr::with_seed(35L, {
  ns <- 50
  prep <- structure(list(), class = "brmsprep")
  prep$dpars <- list(sigma = rchisq(ns, 3), beta = rchisq(ns, 3),
                     mu = matrix(rnorm(ns * 2), ncol = 2),
                     alpha = rnorm(ns), ndt = 1)
  prep$data <- list(Y = rlnorm(ns))
  y1 <- prep$data$Y[1]
  cat("Y[1] =", num(y1), "(shifted_lognormal is finite only above ndt = 1)\n")
  d <- prep$dpars
  cmp("lognormal (log_lik:50)", brms_adapt_log_lik("lognormal", 1, prep),
      dlnorm(y1, d$mu[, 1], d$sigma, log = TRUE))
  cmp("shifted_lognormal (:57)",
      brms_adapt_log_lik("shifted_lognormal", 1, prep),
      brms::dshifted_lnorm(y1, d$mu[, 1], d$sigma, d$ndt, log = TRUE))
  cmp("exgaussian (log_lik:64)", brms_adapt_log_lik("exgaussian", 1, prep),
      brms::dexgaussian(y1, d$mu[, 1], d$sigma, d$beta, log = TRUE))
  cmp("skew_normal (log_lik:71)", brms_adapt_log_lik("skew_normal", 1, prep),
      brms::dskew_normal(y1, d$mu[, 1], d$sigma, d$alpha, log = TRUE))
})

cat("\n## 2. Block 'count and survival', seed 233 (the generated file's)\n")
withr::with_seed(233L, {
  ns <- 25
  nobs <- 10
  trials <- sample(10:30, nobs, replace = TRUE)
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(eta = matrix(rnorm(ns * nobs), ncol = nobs),
                     shape = rgamma(ns, 4), xi = runif(ns, -1, 0.5),
                     phi = rgamma(ns, 1))
  prep$dpars$sigma <- 1 / prep$dpars$shape
  prep$dpars$nu <- prep$dpars$shape + 1
  prep$data <- list(Y = rbinom(nobs, size = trials,
                               prob = rbeta(nobs, 1, 1)), trials = trials)
  i <- sample(nobs, 1)
  cat("i =", i, " Y[i] =", prep$data$Y[i], " trials[i] =", trials[i], "\n")
  y <- prep$data$Y[i]
  d <- prep$dpars
  prep$dpars$mu <- inv_logit(d$eta)
  mu <- prep$dpars$mu[, i]
  cmp("binomial (log_lik:260)", brms_adapt_log_lik("binomial", i, prep),
      dbinom(y, trials[i], mu, log = TRUE))
  cat("extraDistr installed:", requireNamespace("extraDistr", quietly = TRUE),
      "\n")
  r <- tryCatch(brms::dbeta_binomial(y, trials[i], mu, d$phi, log = TRUE),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("brms::dbeta_binomial (the log_lik:262 reference):", r[1], "\n")
  # the beta-binomial in closed form, as an independent check only; the
  # ledger row keeps brms's reference and stays unmeasured
  a <- mu * d$phi
  b <- (1 - mu) * d$phi
  cf <- lchoose(trials[i], y) + lbeta(y + a, trials[i] - y + b) - lbeta(a, b)
  cmp("beta_binomial vs closed form",
      brms_adapt_log_lik("beta_binomial", i, prep), cf)
  prep$dpars$mu <- exp(d$eta)
  mu <- prep$dpars$mu[, i]
  cmp("poisson (log_lik:278)", brms_adapt_log_lik("poisson", i, prep),
      dpois(y, mu, log = TRUE))
  cmp("negbinomial (log_lik:285)", brms_adapt_log_lik("negbinomial", i, prep),
      dnbinom(y, mu = mu, size = d$shape, log = TRUE))
  cmp("geometric (log_lik:295)", brms_adapt_log_lik("geometric", i, prep),
      dnbinom(y, mu = mu, size = 1, log = TRUE))
  # brms's own reference is internal (brms:::dcom_poisson), so the ledger
  # row is not run; read here for information only
  cmp("com_poisson vs brms:::", brms_adapt_log_lik("com_poisson", i, prep),
      brms:::dcom_poisson(y, mu, d$shape, log = TRUE))
  cmp("exponential (log_lik:308)", brms_adapt_log_lik("exponential", i, prep),
      dexp(y, 1 / mu, log = TRUE))
  cmp("gamma (log_lik:316)", brms_adapt_log_lik("gamma", i, prep),
      dgamma(y, shape = d$shape, scale = mu / d$shape, log = TRUE))
  cmp("inverse.gaussian (:331)",
      brms_adapt_log_lik("inverse.gaussian", i, prep),
      brms::dinv_gaussian(y, mu, d$shape, log = TRUE))
  cmp("weibull (log_lik:339)", brms_adapt_log_lik("weibull", i, prep),
      dweibull(y, d$shape, mu / gamma(1 + 1 / d$shape), log = TRUE))
  for (f in c("discrete_weibull", "negbinomial2", "frechet",
              "gen_extreme_value")) {
    cat(f, ": ", tryCatch(brmsfamily(f)$family, error = function(e) {
      substr(conditionMessage(e), 1, 45)
    }), "\n", sep = "")
  }
})

cat("\n## 3. Block 'bernoulli and beta', seed 352\n")
withr::with_seed(352L, {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(mu = inv_logit(matrix(rnorm(ns * nobs * 2),
                                           ncol = nobs * 2)),
                     phi = rgamma(ns, 4))
  prep$data <- list(Y = sample(0:1, nobs, replace = TRUE))
  i <- sample(1:nobs, 1)
  mu <- prep$dpars$mu[, i]
  cmp("bernoulli (log_lik:367)", brms_adapt_log_lik("bernoulli", i, prep),
      dbinom(prep$data$Y[i], 1, mu, log = TRUE))
  prep$data <- list(Y = rbeta(nobs, 1, 1))
  cmp("beta (log_lik:375)", brms_adapt_log_lik("beta", i, prep),
      dbeta(prep$data$Y[i], mu * prep$dpars$phi, (1 - mu) * prep$dpars$phi,
            log = TRUE))
})

cat("\n## 4. log_lik:392, von_mises with cens, seed 378\n")
withr::with_seed(378L, {
  ns <- 15
  nobs <- 10
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(mu = 2 * atan(matrix(rnorm(ns * nobs * 2),
                                          ncol = nobs * 2)),
                     kappa = rgamma(ns, 4))
  prep$data <- list(Y = runif(nobs, -pi, pi))
  i <- sample(seq_len(nobs), 1)
  cmp("von_mises vs brms::dvon_mises",
      brms_adapt_log_lik("von_mises", i, prep),
      brms::dvon_mises(prep$data$Y[i], prep$dpars$mu[, i], prep$dpars$kappa,
                       log = TRUE))
  prep$data$cens <- sample(-1:1, nobs, TRUE)
  cat("i =", i, " cens[i] =", prep$data$cens[i], "\n")
  cat("adapter:", tryCatch(length(brms_adapt_log_lik("von_mises", i, prep)),
                           error = function(e) conditionMessage(e)), "\n")
  ll_b <- brms:::log_lik_von_mises(i, prep = prep)
  cat("brms:::log_lik_von_mises: length", length(ll_b), " first",
      num(ll_b[1]), "\n")
  vm <- von_mises()
  cat("frmtmb von_mises: lcdf", is.function(vm$lcdf), " lccdf",
      is.function(vm$lccdf), " accepts_aterms:",
      paste(vm$accepts_aterms, collapse = ", "), "\n")
  dd <- data.frame(y = prep$data$Y, c = prep$data$cens)
  r <- tryCatch(frm(bf(y | cens(c) ~ 1), family = von_mises(), data = dd),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("frm(y | cens(c) ~ 1, von_mises()):\n  ",
      if (is.character(r)) r else "fits", "\n")
})

cat("\n## 5. posterior_predict truncation reaches frmtmb's rejection step\n")
withr::with_seed(384L, {
  ns <- 30
  nobs <- 15
  prep <- structure(list(ndraws = ns, nobs = nobs), class = "brmsprep")
  prep$dpars <- list(mu = matrix(rnorm(ns * nobs), ncol = nobs),
                     sigma = rchisq(ns, 3))
  prep$refcat <- 1
  prep$data <- list(lb = sample(-(4:7), nobs, TRUE))
  pr <- sapply(1:nobs, brms_adapt_posterior_predict, fun = "gaussian",
               prep = prep)
  cat("gaussian lb: draws below their row's lb:",
      sum(sweep(pr, 2, prep$data$lb) < 0), "of", length(pr), "\n")
  prep$dpars$mu <- exp(prep$dpars$mu)
  prep$data <- list(ub = sample(70:80, nobs, TRUE))
  pr <- sapply(1:nobs, brms_adapt_posterior_predict, fun = "poisson",
               prep = prep)
  cat("poisson ub: draws above their row's ub:",
      sum(sweep(pr, 2, prep$data$ub) > 0), "of", length(pr),
      "; integer:", all(pr == round(pr)), "\n")
})

cat("\n## 6. com_poisson: brms's mu against frmtmb's compois mu\n")
# brms 2.23.0 dcom_poisson(): p(y) proportional to (mu^y / y!)^shape, so
# its mu is lambda^(1 / shape); frmtmb's compois is RTMB::dcompois2(),
# whose mu is the distribution's mean. Same numbers, two distributions.
y <- 0:400
for (pr in list(c(3, 2), c(3, 0.5), c(10, 3), c(3, 1))) {
  mu <- pr[1]
  sh <- pr[2]
  pb <- exp(brms:::dcom_poisson(y, mu, sh, log = TRUE))
  pf <- exp(frmtmb::row_lpdf(compois(), y, y,
                             list(mu = rep(mu, length(y)),
                                  nu = rep(sh, length(y))), list(), NULL))
  cat(sprintf(paste("mu=%g shape=%g: brms mass %.12f mean %.6f;",
                    "frmtmb mass %.12f mean %.6f; log p(y=3) brms %.6f",
                    "frmtmb %.6f\n"),
              mu, sh, sum(pb), sum(y * pb), sum(pf), sum(y * pf),
              log(pb[4]), log(pf[4])))
}

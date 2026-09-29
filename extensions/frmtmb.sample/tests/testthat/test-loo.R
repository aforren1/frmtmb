# log_lik(), loo(), waic() and bayes_R2() on frm_sample() draws: the
# statistical core of the brmsfit post-processing surface. The claim
# under test is that a column of log_lik() is the SAME number the
# objective sums, so the tests are against closed-form densities at a
# known parameter vector rather than against a stored value.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

ll_data <- function(seed = 9, n = 60L, ng = 6L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(ng, 0, 0.5)[dd$g], 1)
  dd
}

# one sampled model per file: the draws are read by most tests here and
# resampling for each would dominate the runtime
ll_case <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      dd <- ll_data()
      fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
      ds <- suppressWarnings(suppressMessages(
        frm_sample(fit, chains = 2, iter = 400, refresh = 0, seed = 1)))
      cache <<- list(dd = dd, fit = fit, ds = ds)
    }
    cache
  }
})

# a draws object with no sampler behind it, for the refusals that are
# decided by the model's STRUCTURE alone: assembling the objective is
# cheap and neither fitting nor sampling adds anything to the test
fake_draws <- function(fit, n = 4L) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  structure(list(stanfit = NULL,
                 draws = matrix(0, n, length(lab),
                                dimnames = list(NULL, lab)),
                 fit = fit),
            class = "frmtmb_draws")
}

## ---- log_lik() is the conditional density, exactly ------------------

test_that("log_lik() equals the closed-form gaussian density at each draw", {
  cs <- ll_case()
  ll <- log_lik(cs$ds)
  expect_equal(dim(ll), c(nrow(cs$ds$draws), nrow(cs$dd)))

  idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
  for (i in c(1L, 7L, 200L, nrow(cs$ds$draws))) {
    sh <- frmtmb.sample:::draws_fit_at(cs$ds, i, idx)
    b <- sh$estimates[["b"]]
    beta <- sh$estimates$beta
    # conditional on THIS draw's own group-level values, which is what
    # brms means by log_lik and what tmbstan makes available
    mu <- beta[1] + beta[2] * cs$dd$x + b[as.integer(cs$dd$g)]
    sig <- exp(sh$estimates$betad[1])
    expect_equal(ll[i, ], stats::dnorm(cs$dd$y, mu, sig, log = TRUE),
                 tolerance = 1e-12)
  }
})

test_that("log_lik() row sums are the taped likelihood at the same draw", {
  cs <- ll_case()
  ll <- log_lik(cs$ds)
  idx <- frmtmb.sample:::draws_par_index(cs$ds$fit)
  i <- 42L
  sh <- frmtmb.sample:::draws_fit_at(cs$ds, i, idx)
  # the objective is -(row densities + random-effect prior), so removing
  # the block density leaves exactly the row sum
  pars <- sh$estimates
  nll <- frmtmb::build_objective(sh$frame)(pars)
  bk <- sh$frame$re_blocks[[1L]]
  re_prior <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll(
    pars[["b"]][bk$b_idx], pars$theta[bk$theta_idx], bk)
  expect_equal(sum(ll[i, ]), as.numeric(-nll) - as.numeric(re_prior),
               tolerance = 1e-10)
})

test_that("log_lik() carries weights, truncation and censoring like the objective", {
  skip_sampler()
  set.seed(4)
  dd <- data.frame(x = stats::rnorm(40))
  dd$y <- stats::rnorm(40, 1 + 0.5 * dd$x, 1)
  dd$w <- rep(c(1, 2), 20)

  # weights multiply the row's contribution (brms's log_lik_weight)
  fw <- frm(bf(y | weights(w) ~ x), family = gaussian(), data = dd)
  dw <- suppressWarnings(suppressMessages(
    frm_sample(fw, chains = 1, iter = 200, refresh = 0, seed = 2)))
  llw <- log_lik(dw)
  idx <- frmtmb.sample:::draws_par_index(dw$fit)
  sh <- frmtmb.sample:::draws_fit_at(dw, 5L, idx)
  mu <- sh$estimates$beta[1] + sh$estimates$beta[2] * dd$x
  sg <- exp(sh$estimates$betad[1])
  expect_equal(llw[5L, ], dd$w * stats::dnorm(dd$y, mu, sg, log = TRUE),
               tolerance = 1e-12)

  # truncation renormalizes by the window mass
  dt <- dd[dd$y > -1, ]
  ft <- frm(bf(y | trunc(lb = -1) ~ x), family = gaussian(), data = dt)
  dts <- suppressWarnings(suppressMessages(
    frm_sample(ft, chains = 1, iter = 200, refresh = 0, seed = 2)))
  llt <- log_lik(dts)
  idx <- frmtmb.sample:::draws_par_index(dts$fit)
  sh <- frmtmb.sample:::draws_fit_at(dts, 5L, idx)
  mu <- sh$estimates$beta[1] + sh$estimates$beta[2] * dt$x
  sg <- exp(sh$estimates$betad[1])
  expect_equal(llt[5L, ],
               stats::dnorm(dt$y, mu, sg, log = TRUE) -
                 stats::pnorm(-1, mu, sg, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-10)

  # right censoring replaces the density by the survival probability
  dc <- dd
  dc$cens <- ifelse(dc$y > 1.5, "right", "none")
  dc$y <- pmin(dc$y, 1.5)
  fc <- frm(bf(y | cens(cens) ~ x), family = gaussian(), data = dc)
  dcs <- suppressWarnings(suppressMessages(
    frm_sample(fc, chains = 1, iter = 200, refresh = 0, seed = 2)))
  llc <- log_lik(dcs)
  idx <- frmtmb.sample:::draws_par_index(dcs$fit)
  sh <- frmtmb.sample:::draws_fit_at(dcs, 5L, idx)
  mu <- sh$estimates$beta[1] + sh$estimates$beta[2] * dc$x
  sg <- exp(sh$estimates$betad[1])
  want <- stats::dnorm(dc$y, mu, sg, log = TRUE)
  r <- dc$cens == "right"
  want[r] <- log(stats::pnorm(dc$y[r], mu[r], sg, lower.tail = FALSE))
  expect_equal(llc[5L, ], want, tolerance = 1e-10)
})

test_that("log_lik() takes trials() through the family, and ndraws thins", {
  skip_sampler()
  set.seed(6)
  dd <- data.frame(x = stats::rnorm(40), n = sample(5:12, 40, TRUE))
  dd$s <- stats::rbinom(40, dd$n, stats::plogis(0.2 + 0.6 * dd$x))
  fit <- frm(bf(s | trials(n) ~ x), family = binomial(), data = dd)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 2)))
  ll <- log_lik(ds)
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  sh <- frmtmb.sample:::draws_fit_at(ds, 3L, idx)
  p <- stats::plogis(sh$estimates$beta[1] + sh$estimates$beta[2] * dd$x)
  expect_equal(ll[3L, ], stats::dbinom(dd$s, dd$n, p, log = TRUE),
               tolerance = 1e-10)

  thin <- log_lik(ds, ndraws = 20)
  expect_equal(dim(thin), c(20L, nrow(dd)))
  expect_null(attr(thin, "chain_id"))
})

test_that("an observation-level mixture keeps its per-row column", {
  skip_sampler()
  set.seed(4)
  dd <- data.frame(y = c(stats::rnorm(50, -2), stats::rnorm(50, 3)))
  fit <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()),
             data = dd)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 2)))
  # only a GROUP-level mixture loses the per-observation column; an
  # ordinary mixture density is a per-row log-sum-exp and factors
  ll <- log_lik(ds)
  expect_equal(dim(ll), c(ndraws(ds), 100L))

  # no random effects here, so the objective IS the negative row sum:
  # an end-to-end check of the whole composition with nothing
  # reconstructed by hand
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  for (i in c(6L, 50L)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
    nll <- frmtmb::build_objective(sh$frame)(sh$estimates)
    expect_equal(sum(ll[i, ]), as.numeric(-nll), tolerance = 1e-10)
  }
})

test_that("log_lik() on a rescor model is the joint density per row", {
  skip_sampler()
  skip_if_not_installed("mvtnorm")
  set.seed(8)
  n <- 50
  Z <- mvtnorm::rmvnorm(n, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2))
  dd <- data.frame(x = stats::rnorm(n), y1 = Z[, 1], y2 = Z[, 2])
  dd$y1 <- dd$y1 + 0.4 * dd$x
  fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE),
             family = gaussian(), data = dd)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 2)))
  ll <- log_lik(ds)
  # one column per OBSERVATION, not per response: the row's two
  # responses share one bivariate density
  expect_equal(ncol(ll), n)

  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  sh <- frmtmb.sample:::draws_fit_at(ds, 4L, idx)
  dp <- frmtmb::eval_dpars(sh)
  C <- frmtmb::us_chol_cor(sh$estimates[["thetar"]], 2L)
  s1 <- rep(as.numeric(dp$y1$sigma), length.out = n)
  s2 <- rep(as.numeric(dp$y2$sigma), length.out = n)
  mu <- cbind(rep(as.numeric(dp$y1$mu), length.out = n),
              rep(as.numeric(dp$y2$mu), length.out = n))
  want <- vapply(seq_len(n), function(i) {
    D <- diag(c(s1[i], s2[i]))
    mvtnorm::dmvnorm(c(dd$y1[i], dd$y2[i]), mu[i, ], D %*% C %*% D,
                     log = TRUE)
  }, numeric(1))
  expect_equal(ll[4L, ], want, tolerance = 1e-8)
})

test_that("the pointwise density of student rescor and mv ordinal fits", {
  # lane mv, 2026-09-25: student() rescor is a multivariate t with one
  # shared nu, and a multivariate model's thresholds are namespaced by
  # response. The rows read at the ML point sum to logLik(), which is
  # the taped objective; before, the rescor rows were the gaussian
  # density and an ordinal response read no thresholds at all
  skip_if_not_installed("mvtnorm")
  set.seed(9)
  n <- 80
  E <- mvtnorm::rmvt(n, matrix(c(1, 0.4, 0.4, 1), 2), df = 4)
  dd <- data.frame(x = stats::rnorm(n))
  dd$y1 <- 1 + dd$x + E[, 1]
  dd$y2 <- E[, 2]
  dd$o <- cut(dd$x + stats::rlogis(n), c(-Inf, -1, 0, 1, Inf),
              labels = FALSE)
  fs <- frm(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + student(),
            data = dd)
  ll <- frmtmb.sample:::draws_row_loglik(fs, NULL)
  L <- as.numeric(logLik(fs))
  expect_lt(abs(sum(ll) - L), 1e3 * .Machine$double.eps * abs(L))
  dp <- frmtmb::eval_dpars(fs)
  C <- frmtmb::us_chol_cor(fs$estimates[["thetar"]], 2L)
  D <- diag(c(dp$y1$sigma[1], dp$y2$sigma[1]))
  want <- mvtnorm::dmvt(c(dd$y1[3], dd$y2[3]),
                        delta = c(dp$y1$mu[3], dp$y2$mu[3]),
                        sigma = D %*% C %*% D, df = dp$y2$nu[1], log = TRUE)
  expect_lt(abs(ll[3] - want), 1e3 * .Machine$double.eps * abs(want))
  fo <- frm(bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian(), data = dd)
  lo <- frmtmb.sample:::draws_row_loglik(fo, NULL)
  L <- as.numeric(logLik(fo))
  expect_lt(abs(sum(lo) - L), 1e3 * .Machine$double.eps * abs(L))
})

## ---- loo() and waic() -----------------------------------------------

test_that("loo() and waic() delegate to the loo package with r_eff", {
  cs <- ll_case()
  skip_if_not_installed("loo")
  l <- suppressWarnings(loo(cs$ds))
  expect_s3_class(l, "loo")
  expect_true(all(is.finite(l$estimates[, "Estimate"])))
  # p_loo is an effective parameter count: eleven parameters here, of
  # which six are shrunk group-level ones
  expect_gt(l$estimates["p_loo", "Estimate"], 1)
  expect_lt(l$estimates["p_loo", "Estimate"], 11)

  ll <- log_lik(cs$ds)
  ref <- suppressWarnings(loo::loo.matrix(
    ll, r_eff = loo::relative_eff(exp(ll),
                                  chain_id = frmtmb.sample:::draws_chain_id(cs$ds))))
  expect_equal(l$estimates, ref$estimates)

  w <- suppressWarnings(waic(cs$ds))
  expect_s3_class(w, "waic")
  # elpd_waic and elpd_loo answer the same question and differ by far
  # less than their own standard error on a well-behaved model
  expect_lt(abs(w$estimates["elpd_waic", 1] - l$estimates["elpd_loo", 1]),
            l$estimates["elpd_loo", 2])

  # thinning drops the chain structure r_eff needs, and says so by
  # giving the matrix no chain_id
  expect_null(attr(log_lik(cs$ds, ndraws = 50), "chain_id"))
  expect_s3_class(suppressWarnings(loo(cs$ds, ndraws = 50)), "loo")
})

test_that("loo_compare() ranks the model that generated the data first", {
  skip_sampler()
  skip_if_not_installed("loo")
  dd <- ll_data(seed = 21, n = 80L, ng = 8L)
  dd$y <- dd$y + 1.5 * dd$x        # a strong x effect to detect
  s1 <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
               chains = 1, iter = 600, refresh = 0, seed = 5)))
  s0 <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ 1 + (1 | g)), family = gaussian(), data = dd,
               chains = 1, iter = 600, refresh = 0, seed = 5)))
  full <- suppressWarnings(loo(s1))
  null <- suppressWarnings(loo(s0))

  # loo reports the model names itself, and where it puts them has moved
  # across versions (a `model` column in current loo, row names before)
  cmp_names <- function(cmp) {
    if (!is.null(cmp$model)) as.character(cmp$model) else rownames(cmp)
  }

  # handed draws, loo_compare() computes the criterion itself and names
  # the models after the arguments
  cmp <- suppressWarnings(loo_compare(s1, s0))
  expect_s3_class(cmp, "compare.loo")
  expect_equal(cmp_names(cmp)[1L], "s1")
  expect_equal(unname(cmp[1L, "elpd_diff"]), 0)
  expect_lt(cmp[2L, "elpd_diff"], -2 * cmp[2L, "se_diff"])
  expect_gt(full$estimates["elpd_loo", 1L],
            null$estimates["elpd_loo", 1L])

  # handed criteria, it is loo's own comparison: frmtmb's generic must
  # not mask it for a script that computed loo() first
  expect_equal(suppressWarnings(loo_compare(full, null)),
               suppressWarnings(loo::loo_compare(full, null)))
  expect_equal(cmp_names(suppressWarnings(
    loo_compare(s1, s0, model_names = c("with_x", "without_x"))))[1L],
    "with_x")
  expect_error(loo_compare(s1, full), "every argument has to be")

  # waic is available on the same path
  expect_s3_class(suppressWarnings(loo_compare(s1, s0,
                                               criterion = "waic")),
                  "compare.loo")
})

test_that("psis() smooths the leave-one-out importance ratios", {
  cs <- ll_case()
  skip_if_not_installed("loo")
  p <- suppressWarnings(psis(cs$ds))
  expect_s3_class(p, "psis")
  ll <- log_lik(cs$ds)
  ref <- suppressWarnings(loo::psis(-ll, r_eff = loo::relative_eff(
    exp(ll), chain_id = frmtmb.sample:::draws_chain_id(cs$ds))))
  expect_equal(loo::pareto_k_values(p), loo::pareto_k_values(ref))
})

test_that("the deprecated capitalized spellings name their replacement", {
  cs <- ll_case()
  expect_error(LOO(cs$ds), "deprecated brms spelling")
  expect_error(WAIC(cs$ds), "Use waic\\(x\\)")
})

## ---- bayes_R2() -----------------------------------------------------

test_that("bayes_R2() is the residual-based estimator of Gelman et al", {
  cs <- ll_case()
  R2 <- bayes_R2(cs$ds, summary = FALSE)
  expect_equal(dim(R2), c(nrow(cs$ds$draws), 1L))
  expect_true(all(R2 > 0 & R2 < 1))

  ep <- posterior_epred(cs$ds)
  vp <- apply(ep, 1, stats::var)
  ve <- apply(sweep(ep, 2, cs$dd$y), 1, stats::var)
  expect_equal(as.numeric(R2), vp / (vp + ve), tolerance = 1e-12)

  s <- bayes_R2(cs$ds)
  expect_equal(colnames(s), c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_equal(rownames(s), "R2")
  expect_equal(s[1, "Estimate"], mean(R2), tolerance = 1e-12)
})

test_that("bayes_R2() refuses an outcome with no residual variance", {
  skip_sampler()
  set.seed(3)
  dd <- data.frame(x = stats::rnorm(60))
  dd$o <- factor(cut(dd$x + stats::rnorm(60), 3), ordered = TRUE)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(o ~ x), family = cumulative(), data = dd,
               chains = 1, iter = 200, refresh = 0, seed = 2)))
  expect_error(bayes_R2(ds), "ordinal or categorical")
})

## ---- refusals -------------------------------------------------------

test_that("log_lik() refuses likelihoods with no per-observation column", {
  # R-side residual correlation: the density is a joint one per group
  set.seed(2)
  da <- data.frame(subj = factor(rep(1:8, each = 5)),
                   week = rep(1:5, 8))
  da$x <- stats::rnorm(40)
  da$y <- stats::rnorm(40, da$x)
  ua <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)), family = gaussian(),
            data = da, dry_run = "objective")
  expect_error(log_lik(fake_draws(ua)), "residual correlation")
  expect_error(loo(fake_draws(ua)), "residual correlation")

  # a group-level mixture: the smallest unit is the group
  set.seed(5)
  dm <- data.frame(id = factor(rep(1:20, each = 4)))
  dm$y <- stats::rnorm(80, rep(c(0, 4), each = 40))
  um <- frm(bf(y ~ 1),
            family = mixture(gaussian(), gaussian(), groups = ~id),
            data = dm, dry_run = "objective")
  expect_error(log_lik(fake_draws(um)), "group-level mixture")
})

test_that("log_lik() refuses in-model imputation and unknown responses", {
  set.seed(7)
  dd <- data.frame(z = stats::rnorm(50))
  dd$x <- dd$z + stats::rnorm(50)
  dd$y <- 1 + dd$x + stats::rnorm(50)
  dd$x[c(3, 8, 20)] <- NA
  um <- frm(bf(y ~ mi(x) + z) + gaussian() + bf(x | mi() ~ z) + gaussian(),
            data = dd, dry_run = "objective")
  expect_error(log_lik(fake_draws(um)), "in-model imputation")

  uf <- frm(bf(y ~ z), family = gaussian(), data = dd[!is.na(dd$x), ],
            dry_run = "objective")
  expect_error(log_lik(fake_draws(uf), resp = "nope"),
               "names no response")
})

test_that("log_lik() refuses laplace-marginalized draws", {
  dd <- ll_data()
  uf <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
            dry_run = "objective")
  lab <- c(frmtmb::brms_par_labels(uf, include_random = FALSE), "lp__")
  ds <- structure(list(stanfit = NULL,
                       draws = matrix(0, 4L, length(lab),
                                      dimnames = list(NULL, lab)),
                       fit = uf),
                  class = "frmtmb_draws")
  expect_error(log_lik(ds), "integrates them out")
})

test_that("the refit and marginal-likelihood methods refuse with a reason", {
  dd <- ll_data()
  uf <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
            dry_run = "objective")
  fd <- fake_draws(uf)
  expect_error(loo_moment_match(fd), "moment-matches")
  expect_error(loo_subsample(fd), "cheap approximation")
  expect_error(reloo(fd), "once per observation")
  expect_error(kfold(fd), "K refits")
  expect_error(bridge_sampler(fd), "integral of the likelihood")
  expect_error(bayes_factor(fd), "ratio of the marginal likelihoods")
  expect_error(post_prob(fd), "posterior model probability")
})

test_that("the draws-only estimators refuse a fit with a pointer", {
  dd <- ll_data()
  uf <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
            dry_run = "objective")
  class(uf) <- "frmtmb_fit"
  expect_error(loo(uf), "posterior quantity")
  expect_error(waic(uf), "no importance|averages the likelihood")
  expect_error(bayes_R2(uf), "per posterior draw")
  expect_error(LOO(uf), "deprecated brms spelling, and on a")
  expect_error(WAIC(uf), "deprecated brms spelling, and it averages")
  expect_error(expose_functions(uf), "no Stan program to read")
})

test_that("multi-model loo() and the draws-surface odd ends refuse", {
  dd <- ll_data()
  uf <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
            dry_run = "objective")
  fd <- fake_draws(uf)
  expect_error(loo(fd, fd), "one model here")
  expect_error(waic(fd, fd), "one model here")
  expect_error(plot(fd), "mcmc_plot")
  expect_error(update(fd), "no formula to revise")
  expect_error(rescor_matrix(fd), "fitted point estimate")
})

## ---- brms's cov = FALSE ARMA has a column per observation ------------

# brms 2.23.0's model block, transliterated, and its J_lag. Copied from
# frmtmb's tests/testthat/test-autocor-cond.R rather than reached for:
# the claim here is that log_lik() reports brms's one-step density, so
# the reference must not be core's own recursion.
#
# for (n in 1:N) {
#   mu[n] += Err[n, 1:Kma] * ma;
#   err[n] = Y[n] - mu[n];
#   for (i in 1:J_lag[n]) Err[n + 1, i] = err[n + 1 - i];
#   mu[n] += Err[n, 1:Kar] * ar;
# }
arma_brms_mu <- function(mu, Y, J_lag, ar, ma) {
  N <- length(Y)
  Err <- matrix(0, N + 1, max(length(ar), length(ma)))
  err <- numeric(N)
  for (n in seq_len(N)) {
    if (length(ma)) mu[n] <- mu[n] + sum(Err[n, seq_along(ma)] * ma)
    err[n] <- Y[n] - mu[n]
    for (i in seq_len(J_lag[n])) Err[n + 1, i] <- err[n + 1 - i]
    if (length(ar)) mu[n] <- mu[n] + sum(Err[n, seq_along(ar)] * ar)
  }
  mu
}

arma_brms_jlag <- function(g, max_lag) {
  N <- length(g)
  J <- integer(N)
  for (n in seq_len(N - 1L)) {
    ind <- n:max(1L, n + 1L - max_lag)
    J[n] <- sum(g[ind] == g[n + 1L])
  }
  J
}

# ragged groups with interior gaps, rows shuffled: brms counts lags in
# ROWS of the (gr, time) order, so neither the gaps nor the data order
# may reach the answer
arma_data <- function(seed = 21L, ng = 8L, nt = 8L) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(nt), g = factor(seq_len(ng)))
  d <- d[-c(3L, 11L, 12L, 30L, 44L), ]
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(nrow(d))
  d <- d[sample(nrow(d)), ]
  rownames(d) <- NULL
  d
}

# A draws object at CHOSEN parameter values and no sampler behind it.
# The stored matrix carries brms's scale where brms has the same
# parameter (`b_*`, `sigma`, `nu`, `r_*`) and frmtmb's internal one
# where it does not: `theta_1` is a log standard deviation, `thetaac_k`
# the raw ARMA coefficients with AR before MA. A vector of length
# `n` per name gives each draw its own value, which is what makes the
# per-draw loop testable without NUTS.
arma_draws <- function(fit, vals, n = 3L) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  m <- matrix(0, n, length(lab), dimnames = list(NULL, lab))
  for (nm in names(vals)) m[, nm] <- vals[[nm]]
  structure(list(stanfit = NULL, draws = m, fit = fit),
            class = "frmtmb_draws")
}

test_that("log_lik() is brms's one-step density for a cov = FALSE ar()", {
  d <- arma_data()
  fit <- frm(bf(y ~ x + ar(t, g)), family = gaussian(), data = d,
             dry_run = "objective")
  # the third draw sets ar = 0, where the shift must vanish entirely
  ars <- c(0.55, -0.3, 0)
  ds <- arma_draws(fit, list(b_Intercept = 1.2, b_x = 0.4, sigma = 0.9,
                             thetaac_1 = ars))
  ll <- log_lik(ds)
  expect_equal(dim(ll), c(3L, nrow(d)))

  ord <- order(d$g, d$t)
  mu0 <- 1.2 + 0.4 * d$x
  J <- arma_brms_jlag(as.integer(d$g)[ord], 1L)
  for (k in seq_along(ars)) {
    mu <- arma_brms_mu(mu0[ord], d$y[ord], J, ars[k], numeric(0))
    ref <- numeric(nrow(d))
    ref[ord] <- stats::dnorm(d$y[ord], mu, 0.9, log = TRUE)
    expect_equal(ll[k, ], ref, tolerance = 1e3 * .Machine$double.eps)
  }
  # ar = 0 is the plain density, which pins the shift's direction: the
  # draw with ar = 0.55 must NOT be that, by far more than the residual
  plain <- stats::dnorm(d$y, mu0, 0.9, log = TRUE)
  expect_equal(ll[3L, ], plain, tolerance = 1e3 * .Machine$double.eps)
  expect_gt(abs(sum(ll[1L, ]) - sum(plain)), 0.05 * abs(sum(plain)))
})

test_that("a cov = FALSE row sum is the objective at the same draw", {
  d <- arma_data()
  # no random effect: build_objective() is then the row product alone,
  # so the identity needs no correction term
  fit <- frm(bf(y ~ x + arma(t, g, p = 1, q = 1)), family = gaussian(),
             data = d, dry_run = "objective")
  ds <- arma_draws(fit, list(b_Intercept = 1.1, b_x = 0.45,
                             sigma = 0.85, thetaac_1 = c(0.5, -0.2, 0.3),
                             thetaac_2 = c(0.25, 0.4, -0.35)))
  ll <- log_lik(ds)
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  for (k in seq_len(nrow(ds$draws))) {
    sh <- frmtmb.sample:::draws_fit_at(ds, k, idx)
    nll <- frmtmb::build_objective(sh$frame)(sh$estimates)
    expect_equal(sum(ll[k, ]), -as.numeric(nll),
                 tolerance = 1e3 * .Machine$double.eps)
  }
  # the three draws are three different likelihoods, so the identity is
  # not being satisfied by a constant
  expect_gt(stats::sd(rowSums(ll)), 0)
})

test_that("log_lik() carries student, sigma ~ x and (1 | g) under arma()", {
  d <- arma_data()
  fit <- frm(bf(y ~ x + arma(t, g, p = 1, q = 1) + (1 | g)),
             family = student(), data = d, dry_run = "objective")
  rn <- grep("^r_", frmtmb::brms_par_labels(fit), value = TRUE)
  bvals <- seq(-0.5, 0.5, length.out = length(rn))
  v <- as.list(stats::setNames(bvals, rn))
  v$b_Intercept <- 1.1
  v$b_x <- 0.45
  v$sigma <- 0.8
  v$nu <- 7
  v$theta_1 <- log(0.6)
  v$thetaac_1 <- 0.5
  v$thetaac_2 <- 0.25
  ds <- arma_draws(fit, v, n = 2L)
  ll <- log_lik(ds)

  ord <- order(d$g, d$t)
  mu0 <- 1.1 + 0.45 * d$x + bvals[as.integer(d$g)]
  mu <- arma_brms_mu(mu0[ord], d$y[ord],
                     arma_brms_jlag(as.integer(d$g)[ord], 1L), 0.5, 0.25)
  ref <- numeric(nrow(d))
  # brms's student_t_lpdf(nu, mu, sigma), the scaled t
  ref[ord] <- stats::dt((d$y[ord] - mu) / 0.8, df = 7, log = TRUE) -
    log(0.8)
  expect_equal(ll[1L, ], ref, tolerance = 1e3 * .Machine$double.eps)

  # the row product against the objective needs the block's own prior
  # back out: build_objective() carries it and log_lik() must not
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  sh <- frmtmb.sample:::draws_fit_at(ds, 1L, idx)
  bk <- sh$frame$re_blocks[[1L]]
  rp <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll(
    sh$estimates[["b"]][bk$b_idx], sh$estimates[["theta"]][bk$theta_idx],
    bk)
  nll <- frmtmb::build_objective(sh$frame)(sh$estimates)
  expect_equal(sum(ll[1L, ]), -as.numeric(nll) - as.numeric(rp),
               tolerance = 1e3 * .Machine$double.eps)

  # a distributional sigma travels the same way
  fs <- frm(bf(y ~ x + ar(t, g), sigma ~ x), family = gaussian(),
            data = d, dry_run = "objective")
  dss <- arma_draws(fs, list(b_Intercept = 1.2, b_x = 0.4,
                             b_sigma_Intercept = -0.2, b_sigma_x = 0.3,
                             thetaac_1 = 0.55), n = 1L)
  lls <- log_lik(dss)
  sg <- exp(-0.2 + 0.3 * d$x)
  mus <- arma_brms_mu((1.2 + 0.4 * d$x)[ord], d$y[ord],
                      arma_brms_jlag(as.integer(d$g)[ord], 1L), 0.55,
                      numeric(0))
  refs <- numeric(nrow(d))
  refs[ord] <- stats::dnorm(d$y[ord], mus, sg[ord], log = TRUE)
  expect_equal(lls[1L, ], refs, tolerance = 1e3 * .Machine$double.eps)
})

test_that("loo() runs on a cov = FALSE fit and cov = TRUE still refuses", {
  skip_if_not_installed("loo")
  d <- arma_data()
  fit <- frm(bf(y ~ x + ar(t, g)), family = gaussian(), data = d,
             dry_run = "objective")
  # a spread of ar values so PSIS has something to smooth
  set.seed(5)
  n <- 40L
  ds <- arma_draws(fit, list(
    b_Intercept = stats::rnorm(n, 1.2, 0.1),
    b_x = stats::rnorm(n, 0.4, 0.1),
    sigma = exp(stats::rnorm(n, log(0.9), 0.05)),
    thetaac_1 = stats::rnorm(n, 0.5, 0.1)), n = n)
  lo <- suppressWarnings(loo(ds))
  expect_true(is.finite(lo$estimates["elpd_loo", "Estimate"]))
  expect_true(is.finite(lo$estimates["p_loo", "Estimate"]))
  expect_equal(length(lo$pointwise[, "elpd_loo"]), nrow(d))
  expect_true(all(is.finite(lo$pointwise[, "elpd_loo"])))
  # a per-observation matrix carries no unit, so loo() says nothing
  expect_no_message(suppressWarnings(loo(ds)))
  wa <- suppressWarnings(waic(ds))
  expect_true(is.finite(wa$estimates["elpd_waic", "Estimate"]))

  # the covariance form keeps its refusal, and the message names the
  # form that does factor
  uc <- frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = gaussian(),
            data = d, dry_run = "objective")
  expect_error(log_lik(fake_draws(uc)), "residual correlation MATRIX")
  expect_error(log_lik(fake_draws(uc)), "cov = FALSE form of ar()",
               fixed = TRUE)
  expect_error(loo(fake_draws(uc)), "residual correlation MATRIX")
  # cosy() and unstr() have no cov = FALSE form, so the message must not
  # offer one
  uu <- frm(bf(y ~ x + unstr(t, g)), family = gaussian(), data = d,
            dry_run = "objective")
  expect_error(log_lik(fake_draws(uu)), "residual correlation MATRIX")
  expect_error(log_lik(fake_draws(uu)), "unstr(t, g)", fixed = TRUE)
  msg <- conditionMessage(tryCatch(log_lik(fake_draws(uu)),
                                   error = identity))
  expect_false(grepl("cov = FALSE", msg, fixed = TRUE))

  # mi() on the response stays refused for its own reason, which is not
  # about autocorrelation
  d2 <- d
  d2$y[c(4L, 9L)] <- NA
  um <- frm(bf(y | mi() ~ x + ar(t, g)) + gaussian(), data = d2,
            dry_run = "objective")
  expect_error(log_lik(fake_draws(um)), "in-model imputation")
})

test_that("log_lik() matches brms row by row on a cov = FALSE ARMA", {
  skip_unless_brms_fit()
  # no skip_sampler(): the draws are a fixed matrix with stanfit = NULL,
  # so tmbstan is not on this block's path at all
  d <- arma_data(seed = 22L, ng = 6L, nt = 8L)
  bfit <- suppressWarnings(suppressMessages(
    brms::brm(brms::bf(y ~ x + arma(t, g, p = 1, q = 1)), data = d,
              family = brms::brmsfamily("gaussian"), chains = 1,
              iter = 300, warmup = 200, refresh = 0, seed = 5,
              silent = 2)))
  fit <- frm(bf(y ~ x + arma(t, g, p = 1, q = 1)), family = gaussian(),
             data = d, dry_run = "objective")
  # brms's own draws, written into a frmtmb draws object by name: the
  # two packages then report the same quantity at the same parameter
  # vector, and nothing rests on the two samplers agreeing
  bm <- as.matrix(brms::as_draws_matrix(bfit))
  ds <- arma_draws(fit, list(b_Intercept = bm[, "b_Intercept"],
                             b_x = bm[, "b_x"],
                             sigma = bm[, "sigma"],
                             thetaac_1 = bm[, "ar[1]"],
                             thetaac_2 = bm[, "ma[1]"]),
                   n = nrow(bm))
  # the transplant is checked before the claim rests on it: the one-step
  # mean reads every transplanted parameter
  expect_lt(max(abs(posterior_epred(ds) - brms::posterior_epred(bfit))),
            1e-10 * stats::sd(d$y))
  llb <- brms::log_lik(bfit)
  llf <- log_lik(ds)
  expect_equal(dim(llf), dim(llb))
  # brms's log_lik undoes its own order(gr, time) sort (reorder_obs), so
  # both matrices are in the user's row order
  expect_lt(max(abs(llf - llb)), 1e-10 * stats::sd(llb))
})

## ---- brms cross-check (opt-in; compiles Stan) ------------------------

test_that("loo() and log_lik() agree with brms on the same model", {
  skip_unless_brms_fit()
  skip_if_not_installed("loo")
  skip_sampler()
  set.seed(11)
  n <- 80L; ng <- 8L
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x +
                         stats::rnorm(ng, 0, 0.6)[dd$g], 1)

  # the formula route applies brms's own default priors, so the two
  # posteriors are the same target and the comparison is of estimators,
  # not of models
  ds <- suppressWarnings(suppressMessages(
    frm_sample(bf(y ~ x + (1 | g)), family = gaussian(), data = dd,
               chains = 2, iter = 2000, refresh = 0, seed = 3)))
  bfit <- suppressMessages(brms::brm(y ~ x + (1 | g), data = dd,
                                     family = brms::brmsfamily("gaussian"),
                                     chains = 2, iter = 2000, refresh = 0,
                                     seed = 3, silent = 2))

  fl <- suppressWarnings(loo(ds))
  bl <- suppressWarnings(brms::loo(bfit))
  # TOLERANCE: both elpds are Monte Carlo estimates from independent
  # chains of the same posterior, so the honest yardstick is their own
  # reported standard error, not a fixed number of digits. Two
  # independent estimates agree when their difference is small against
  # the pooled SE; 0.5 * pooled SE is a deliberately tight version of
  # that (the SE here is ~7 elpd units on 80 observations).
  se <- sqrt(fl$estimates["elpd_loo", 2]^2 + bl$estimates["elpd_loo", 2]^2)
  expect_lt(abs(fl$estimates["elpd_loo", 1] - bl$estimates["elpd_loo", 1]),
            0.5 * se)
  # p_loo is a much better determined quantity and agrees far more
  # closely: same effective parameter count to within a parameter
  expect_lt(abs(fl$estimates["p_loo", 1] - bl$estimates["p_loo", 1]), 1)

  # the log_lik matrix itself: column means are per-observation expected
  # log-densities, whose MC error is sd/sqrt(ndraws) per column
  cf <- colMeans(log_lik(ds))
  cb <- colMeans(brms::log_lik(bfit))
  expect_lt(max(abs(cf - cb)), 0.1)
  expect_gt(stats::cor(cf, cb), 0.999)

  expect_lt(abs(bayes_R2(ds)[1, "Estimate"] -
                  brms::bayes_R2(bfit)[1, "Estimate"]), 0.05)
})

test_that("a group-unit matrix says so, because loo() cannot", {
  # loo::loo.matrix() prints "Computed from N by K log-likelihood
  # matrix" whatever the columns are, and never sees the attribute
  # log_lik() attaches, so 8 subjects read exactly like 8 observations.
  # The message is the only place a caller learns which they asked for.
  ll <- matrix(rnorm(30), 5, 6)
  attr(ll, "unit") <- "one subject's trial sequence"
  # The METHOD is mocked, not the generic. `log_lik` in this namespace is
  # an active binding to rstantools' generic (R/generic-owners.R), and
  # an assignment to an active binding calls its function with the
  # value, so mocking the generic fails with "unused argument".
  #
  # Mocking the method works here, and the reason bounds when it works.
  # local_mocked_bindings() patches the namespace binding and THIS
  # package's method table, never the owner's. With rstantools loaded,
  # the exported generic is rstantools' and its table still holds the
  # real method, so only a caller INSIDE this namespace sees the mock:
  # UseMethod() searches the calling environment before any table, and
  # loo_matrix() is such a caller. A mock of this shape does not cover a
  # caller at the prompt.
  # .package named so the mock does not depend on a pkgload dev namespace
  fd <- structure(list(), class = "frmtmb_draws")
  local_mocked_bindings(
    log_lik.frmtmb_draws = function(object, ndraws = NULL, resp = NULL,
                                    ...) ll,
    .package = "frmtmb.sample")
  expect_message(loo_matrix(fd, NULL, NULL, "loo()"),
                 "leave-one-out over the 6 units")
  expect_message(loo_matrix(fd, NULL, NULL, "loo()"),
                 "one subject's trial sequence", fixed = TRUE)
  expect_message(loo_matrix(fd, NULL, NULL, "waic()"),
                 "waic() is leave-one-out", fixed = TRUE)
  # a per-observation matrix carries no unit and says nothing
  plain <- matrix(rnorm(30), 5, 6)
  local_mocked_bindings(
    log_lik.frmtmb_draws = function(object, ndraws = NULL, resp = NULL,
                                    ...) plain,
    .package = "frmtmb.sample")
  expect_no_message(loo_matrix(fd, NULL, NULL, "loo()"))
})

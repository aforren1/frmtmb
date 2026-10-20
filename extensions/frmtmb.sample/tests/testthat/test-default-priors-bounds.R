# The defaults the sampler can apply, and the bounds brms declares with
# them (review dev/reviews/2026-10-14-sampleprior.md, findings 1, 3, 4
# and 7).

skip_on_cran()

dpx_mv_data <- function() {
  set.seed(20261114)
  n <- 60
  x <- stats::rnorm(n)
  data.frame(x = x, y1 = 1 + x + stats::rt(n, 5),
             y2 = 2 - x + stats::rt(n, 5))
}

# the objective frm_sample() samples, rebuilt as it builds it, and its
# log prior as that objective less the bare likelihood's
dpx_logprior <- function(fit, prior) {
  rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, prior))
  f1 <- fit
  f1$prior <- rp$effective
  od <- frmtmb.sample:::prior_augmented_obj(f1, rp$ri$entries)
  f0 <- fit
  f0$prior <- NULL
  of <- frmtmb.sample:::prior_augmented_obj(f0, list())
  function(p) of$fn(p) - od$fn(p)
}

test_that("a Student-t rescor model samples with brms's shared nu", {
  skip_sampler()
  d <- dpx_mv_data()
  fam <- frmtmb::brmsfamily("student")
  # through the first lane build every one of these stopped: the shared
  # nu's default has no resp, and resp_missing_refusal() refused it
  forms <- list(mvbind = bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE),
                mvbf = mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE))
  for (nm in names(forms)) {
    ds <- suppressWarnings(suppressMessages(frm_sample(
      forms[[nm]], data = d, family = fam, chains = 1, iter = 200,
      refresh = 0, seed = 1)))
    ps <- as.data.frame(prior_summary(ds))
    expect_identical(ps$prior[ps$class == "nu"], "gamma(2, 0.1)",
                     label = nm)
    expect_true("rescor__y1__y2" %in% colnames(ds$draws), label = nm)
  }
  # the fit route, and a user's prior on the shared nu, which brms
  # writes with no resp
  fit <- suppressWarnings(frm(forms$mvbind, data = d, family = fam))
  ds <- suppressWarnings(suppressMessages(frm_sample(
    fit, chains = 1, iter = 200, refresh = 0, seed = 1,
    prior = set_prior("gamma(3, 0.2)", class = "nu"))))
  ps <- as.data.frame(prior_summary(ds))
  expect_identical(ps$prior[ps$class == "nu"], "gamma(3, 0.2)")
})

test_that("the shared nu's prior in the objective is brms's gamma", {
  d <- dpx_mv_data()
  fit <- suppressWarnings(frm(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE),
                              data = d,
                              family = frmtmb::brmsfamily("student")))
  for (pr in list(NULL, set_prior("gamma(3, 0.2)", class = "nu"))) {
    lp <- dpx_logprior(fit, pr)
    p <- fit$obj$env$last.par.best
    j <- grep("nu", names(fit$frame$par_template$betad))
    expect_length(j, 1L)
    k <- length(fit$frame$par_template$beta) + j
    a <- if (is.null(pr)) c(2, 0.1) else c(3, 0.2)
    # nu = 1 + exp(eta): the gamma density of nu and the log Jacobian
    dens <- function(eta) {
      stats::dgamma(1 + exp(eta), a[1], a[2], log = TRUE) + eta
    }
    p2 <- p
    p2[k] <- p2[k] + 0.7
    expect_equal(unname(lp(p2) - lp(p)), unname(dens(p2[k]) - dens(p[k])),
                 tolerance = 1e-10)
  }
})

test_that("a log-linked ndt default is bounded by the smallest response", {
  set.seed(20261116)
  d <- data.frame(x = stats::rnorm(60))
  d$y <- 0.4 + exp(stats::rnorm(60, -0.5 + 0.2 * d$x, 0.4))
  uf <- frm(bf(y ~ x), family = shifted_lognormal(), data = d,
            dry_run = "objective")
  s <- Filter(function(s) identical(s$dpar, "ndt"),
              unclass(frmtmb.sample:::default_priors_for(uf)))[[1L]]
  # brms declares real<lower=0, upper=min_Y> ndt; the bound is what
  # Stan's transform keeps a draw inside, where the density alone is a
  # wall (82 divergent transitions of 4000 on these data without it)
  expect_identical(s$ub, min(d$y))
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf,
                                                               NULL))$ri
  j <- grep("ndt", names(ri$upper))
  expect_identical(unname(ri$upper[j]), log(min(d$y)))
})

test_that("an identity-linked nu default is bounded below by 1", {
  set.seed(20261016)
  d <- data.frame(x = stats::rnorm(40))
  d$y <- 1 + d$x + stats::rt(40, 5)
  nu_spec <- function(fam) {
    uf <- frm(bf(y ~ x), family = fam, data = d, dry_run = "objective")
    Filter(function(s) identical(s$dpar, "nu"),
           unclass(frmtmb.sample:::default_priors_for(uf)))[[1L]]
  }
  # brms declares real<lower=1> nu whatever the link
  expect_identical(nu_spec(student(link_nu = "identity"))$lb, 1)
  # logm1 holds nu above 1 already, so no bound is written there
  expect_true(is.na(nu_spec(student())$lb))
})

test_that("the fast sd_ path is checked where theta is extreme", {
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(11)
  n <- 120
  dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:12, 10)))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + stats::rnorm(12, 0, 0.7)[dd$g])
  fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3)))
  nc <- frmtmb.sample:::draws_natural_cols(ds$fit)
  im <- frmtmb.sample:::draws_internal_matrix(ds)
  dv <- nc$derived
  th <- im[, "theta_1"]
  top <- max(th)
  # a fast path wrong only at the largest theta, which neither point of
  # the construction's check visits: the extreme-draw check sends every
  # draw through hyp_env_vals() instead
  good <- dv$fast$vals_at
  dv$fast$vals_at <- function(theta) {
    v <- good(theta)
    if (theta[[1L]] >= top) v * (1 + 1e-6) else v
  }
  slow <- dv
  slow$fast <- NULL
  expect_identical(frmtmb.sample:::draws_cov_values(im, ds$fit, dv),
                   frmtmb.sample:::draws_cov_values(im, ds$fit, slow))
})

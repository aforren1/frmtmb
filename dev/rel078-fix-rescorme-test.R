# A multivariate model with residual correlation and me() terms, the one
# model both changes of frmtmb.sample 0.25.0 reach: brms's me() prior
# classes (meanme, sdme, corme) and the joint predictive draw and
# per-response density of set_rescor(TRUE) draws. Each half was tested
# on a model without the other; this is the cross term
# (dev/round-20261029.md, "The semantic overlaps"; probe
# dev/rel078-overlap.R). Seen to fail on frmtmb.sample 0.24.0, which
# refuses the me() classes and returns the first response only, and on
# each half's build without the other.

skip_on_cran()
withr::local_options(mc.cores = 1, .local_envir = teardown_env())

rm_data <- function() {
  set.seed(29)
  n <- 60
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  e <- matrix(stats::rnorm(2 * n), n) %*% chol(matrix(c(1, .7, .7, 1), 2))
  data.frame(xo = x + stats::rnorm(n, 0, 0.3), sx = 0.3,
             zo = z + stats::rnorm(n, 0, 0.25), sz = 0.25,
             y = 1 + 0.8 * x + 0.5 * e[, 1],
             y2 = -0.5 * x + 0.4 * z + 0.6 * e[, 2])
}

rm_form <- function() {
  mvbf(bf(y ~ me(xo, sx)), bf(y2 ~ me(xo, sx) + me(zo, sz))) +
    set_rescor(TRUE)
}

rm_prior <- function() {
  set_prior("normal(0, 2)", class = "meanme") +
    set_prior("exponential(1)", class = "sdme") +
    set_prior("lkj(2)", class = "corme")
}

rm_quiet <- function(expr) suppressWarnings(suppressMessages(expr))

test_that("the sampling-route table lists the me() and rescor rows as brms", {
  tab <- default_prior(rm_form(), data = rm_data(), route = "sample")
  key <- paste(tab$class, tab$coef, tab$group, tab$resp, sep = "|")
  # brms 2.23.0's default_prior() on this formula (dev/rel078-overlap.R,
  # dev/rel078-log/overlap-1.txt): the me() classes carry no resp, once
  # for the model, and rescor is one row
  for (k in c("meanme|||", "meanme|mexo||", "meanme|mezo||", "sdme|||",
              "sdme|mexo||", "sdme|mezo||", "corme|||", "rescor|||")) {
    expect_true(k %in% key, label = k)
  }
  expect_identical(tab$prior[key == "corme|||"], "lkj(1)")
  expect_identical(tab$prior[key == "rescor|||"], "lkj(1)")
  expect_identical(as.numeric(tab$lb[key == "sdme|||"]), 0)
})

rm_draws <- local({
  cache <- NULL
  function() {
    skip_sampler()
    if (is.null(cache)) {
      cache <<- rm_quiet(frm_sample(
        rm_form(), data = rm_data(), family = gaussian(), prior = rm_prior(),
        sample_prior = "yes", chains = 2, iter = 800, seed = 5,
        refresh = 0))
    }
    cache
  }
})

test_that("me() priors and the joint draw meet on one rescor model", {
  ds <- rm_draws()
  m <- as.matrix(ds)
  # the me() hyperparameters and the residual correlation are sampled
  for (v in c("meanme_mexo", "meanme_mezo", "sdme_mexo", "sdme_mezo",
              "corme__mexo__mezo", "rescor__y__y2")) {
    expect_true(v %in% colnames(m), label = v)
  }
  # their prior draws, under brms's names
  pd <- prior_draws(ds)
  expect_true(all(c("meanme_1", "sdme_1", "Lme_1", "rescor") %in% names(pd)))
  # the predictive draw is draws x rows x responses, with both responses
  pp <- posterior_predict(ds)
  ep <- posterior_epred(ds)
  expect_identical(dim(pp), c(nrow(m), nrow(rm_data()), 2L))
  expect_identical(dimnames(pp)[[3L]], c("y", "y2"))
  # drawn jointly: per draw, the mean over rows of the product of the
  # two standardized predictive residuals has expectation the draw's
  # rescor (a sample correlation would be biased by about
  # rho (1 - rho^2) / (2 n)); its error is measured by the spread over
  # draws
  u <- rowMeans((pp[, , 1] - ep[, , 1]) / m[, "sigma_y"] *
                  (pp[, , 2] - ep[, , 2]) / m[, "sigma_y2"])
  dlt <- u - m[, "rescor__y__y2"]
  z <- mean(dlt) / (stats::sd(dlt) / sqrt(length(dlt)))
  expect_lt(abs(z), 4)
  # an independent draw would center r on 0: the bound has power
  z0 <- mean(m[, "rescor__y__y2"]) / (stats::sd(dlt) / sqrt(length(dlt)))
  expect_gt(z0, 8)
})

test_that("log_lik(resp =) of the me() rescor model is each response's own", {
  ds <- rm_draws()
  m <- as.matrix(ds)
  d <- rm_data()
  ep <- posterior_epred(ds)
  ll1 <- log_lik(ds, resp = "y")
  ref1 <- stats::dnorm(matrix(d$y, nrow(m), nrow(d), byrow = TRUE),
                       ep[, , 1], m[, "sigma_y"], log = TRUE)
  expect_lt(max(abs(ll1 - ref1)) / max(abs(ref1)), 1e-12)
  # both, joint: the bivariate normal at the draw's rescor
  ll <- log_lik(ds)
  ll2 <- log_lik(ds, resp = "y2")
  rho <- m[, "rescor__y__y2"]
  z1 <- (matrix(d$y, nrow(m), nrow(d), byrow = TRUE) - ep[, , 1]) /
    m[, "sigma_y"]
  z2 <- (matrix(d$y2, nrow(m), nrow(d), byrow = TRUE) - ep[, , 2]) /
    m[, "sigma_y2"]
  refj <- -log(2 * pi) - log(m[, "sigma_y"]) - log(m[, "sigma_y2"]) -
    0.5 * log(1 - rho^2) -
    (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * (1 - rho^2))
  expect_lt(max(abs(ll - refj)) / max(abs(refj)), 1e-12)
  # and the joint is not the sum of the margins at this correlation
  expect_gt(max(abs(ll - (ll1 + ll2))), 1e-3 * max(abs(refj)))
})

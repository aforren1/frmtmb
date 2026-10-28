# build_objective()'s `prior_only` field, which frmtmb.sample's
# frm_sample(sample_prior = "only") sets on a copy of the frame: brms's
# `if (!prior_only)` around the likelihood. What brms keeps outside it is
# the group-level densities and me()'s measurement and latent densities
# (brms 2.23.0's stancode(sample_prior = "only"), recorded in
# frmtmb.sample's findings for this seam).

po_data <- function() {
  set.seed(1)
  d <- data.frame(x = rnorm(40), g = factor(rep(1:5, 8)))
  d$y <- 1 + d$x + rnorm(5)[d$g] + rnorm(40)
  d
}

test_that("prior_only leaves out the response density and keeps b's", {
  fit <- frm(bf(y ~ x + (1 | g)), data = po_data())
  fr <- fit$frame
  fr[["prior_only"]] <- TRUE
  po <- build_objective(fr)
  full <- build_objective(fit$frame)
  p <- fit$estimates
  # what is left is the group-level density: b ~ normal(0, exp(theta))
  ref <- -sum(stats::dnorm(p$b, 0, exp(p$theta), log = TRUE))
  expect_equal(po(p), ref, tolerance = 1e-12)
  # it does not read the coefficients or sigma, which only the response
  # density reads
  q <- p
  q$beta <- q$beta + 1
  q$betad <- q$betad - 0.5
  expect_identical(po(q), po(p))
  # the full objective still does, and the fit's own frame is untouched
  expect_false(isTRUE(all.equal(full(q), full(p))))
  expect_null(fit$frame[["prior_only"]])
  # the difference is the response's log-likelihood at the estimates
  d <- po_data()
  mu <- drop(fit$frame$linpreds[[1]]$X %*% p$beta) + p$b[as.integer(d$g)]
  ll <- sum(stats::dnorm(d$y, mu, exp(p$betad), log = TRUE))
  expect_equal(full(p) - po(p), -ll, tolerance = 1e-10)
})

test_that("prior_only keeps me()'s measurement and latent densities", {
  set.seed(2)
  n <- 40
  tx <- rnorm(n)
  d <- data.frame(x = tx + rnorm(n, 0, 0.3), sx = 0.3)
  d$y <- 1 + tx + rnorm(n)
  fit <- frm(bf(y ~ me(x, sx)), data = d)
  fr <- fit$frame
  fr[["prior_only"]] <- TRUE
  po <- build_objective(fr)
  p <- fit$estimates
  q <- p
  q$beta <- q$beta + 1
  # the coefficients feed the response only
  expect_identical(po(q), po(p))
  # the latent values have their two densities, as in brms
  q <- p
  q$miss <- q$miss + 0.1
  expect_false(isTRUE(all.equal(po(q), po(p))))
  mu <- p$meanme[[1]]
  s <- exp(p$logsdme[[1]])
  ref <- -sum(stats::dnorm(d$x, p$miss, 0.3, log = TRUE)) -
    sum(stats::dnorm(p$miss, mu, s, log = TRUE))
  expect_equal(po(p), ref, tolerance = 1e-10)
})

test_that("prior_only keeps a mi(sdy) response's measurement density", {
  # brms 2.23.0, stancode(bf(y | mi(sdy) ~ x), sample_prior = "only"):
  # target += normal_lpdf(Y[Jme] | Yl[Jme], noise[Jme]) sits outside
  # `if (!prior_only)`, so the latent responses have a proper density
  set.seed(3)
  n <- 40
  d <- data.frame(x = rnorm(n), sdy = 0.4)
  d$y <- 1 + d$x + rnorm(n)
  fit <- frm(bf(y | mi(sdy) ~ x), data = d)
  fr <- fit$frame
  fr[["prior_only"]] <- TRUE
  po <- build_objective(fr)
  p <- fit$estimates
  ref <- -sum(stats::dnorm(d$y, p$miss, 0.4, log = TRUE))
  expect_equal(po(p), ref, tolerance = 1e-10)
  q <- p
  q$beta <- q$beta + 1
  expect_identical(po(q), po(p))
})

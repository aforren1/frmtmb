# Non-gaussian REML: what frm(REML = TRUE) computes beyond the
# gaussian case, pinned two ways. Consistency: glmmTMB integrates the
# fixed effects out with the same Laplace trick, so the two must agree
# to optimizer precision. Meaning: the objective IS the
# Laplace-approximated integrated likelihood, which for a fixed theta
# equals the profiled nll plus half the log determinant of the
# fixed-effect information over 2 pi, i.e. the Cox-Reid adjusted
# profile likelihood computed automatically.

skip_on_cran()

rng_data <- function(seed = 71, n = 240L, ng = 12L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  dd$yp <- stats::rpois(n, exp(0.4 + 0.3 * dd$x +
                                 stats::rnorm(ng, 0, 0.5)[dd$g]))
  dd$yb <- stats::rbinom(n, 1L,
                         stats::plogis(0.3 + 0.6 * dd$x +
                                         stats::rnorm(ng, 0, 0.7)[dd$g]))
  dd
}

test_that("non-gaussian REML agrees with glmmTMB's", {
  skip_if_not_installed("glmmTMB")
  dd <- rng_data()

  fp <- frm(bf(yp ~ x + (1 | g)), family = poisson(), data = dd,
            REML = TRUE)
  gp <- glmmTMB::glmmTMB(yp ~ x + (1 | g), family = poisson,
                         data = dd, REML = TRUE)
  expect_lt(abs(as.numeric(logLik(fp)) - as.numeric(logLik(gp))), 1e-5)
  expect_lt(abs(sqrt(VarCorr(fp)[[1]][1, 1]) -
                  attr(glmmTMB::VarCorr(gp)$cond$g, "stddev")[[1]]),
            1e-4)

  fb <- frm(bf(yb ~ x + (1 | g)), family = bernoulli(), data = dd,
            REML = TRUE)
  gb <- glmmTMB::glmmTMB(yb ~ x + (1 | g), family = binomial,
                         data = dd, REML = TRUE)
  expect_lt(abs(as.numeric(logLik(fb)) - as.numeric(logLik(gb))), 1e-5)
})

test_that("the REML objective is the Cox-Reid adjusted profile nll", {
  # small model, theta held at the REML estimate: integrate beta out by
  # hand with the Laplace formula and compare against the fit's own
  # objective. The identity is
  #   nll_REML(theta) = nll(beta_hat(theta), theta)
  #                     + 0.5 * log det(H_beta / (2 pi))
  # where H_beta is the Hessian of the beta-profiled objective at its
  # minimum, which is exactly the Cox-Reid adjustment term
  dd <- rng_data(seed = 72, n = 160L, ng = 8L)
  fit <- frm(bf(yp ~ x + (1 | g)), family = poisson(), data = dd,
             REML = TRUE)
  # the ML objective at fixed theta = the REML estimate, as a function
  # of beta only (b still Laplace-integrated): refit ML with theta
  # mapped to the REML estimate via bounds pinning
  th <- fit$estimates$theta
  ml <- suppressWarnings(
    frm(bf(yp ~ x + (1 | g)), family = poisson(), data = dd,
        lower = c(theta_1 = th[[1]]), upper = c(theta_1 = th[[1]]),
        start = list(theta = unname(th))))
  nll_hat <- -as.numeric(logLik(ml))
  # the beta block of the marginal objective's Hessian at the
  # pinned-theta optimum, by central differences of the AD gradient
  # (obj$he is unavailable with random effects); the joint-Laplace
  # determinant factorizes so that this Schur block is exactly the
  # Cox-Reid adjustment matrix
  obj <- ml$obj
  p0 <- ml$opt$par
  bidx <- which(names(p0) == "beta")
  gr_b <- function(b) {
    p <- p0
    p[bidx] <- b
    obj$gr(p)[bidx]
  }
  b0 <- p0[bidx]
  h <- 1e-5
  H <- matrix(NA_real_, 2, 2)
  for (j in 1:2) {
    e <- c(0, 0)
    e[j] <- h
    H[, j] <- (gr_b(b0 + e) - gr_b(b0 - e)) / (2 * h)
  }
  H <- (H + t(H)) / 2
  cr <- nll_hat + 0.5 * (determinant(H)$modulus[1] - 2 * log(2 * pi))
  # the slack is the pinned refit's optimizer tolerance plus finite
  # differences over the inner Laplace solve, measured at ~4e-3 here;
  # a wrong identity (the wrong determinant, a missing 2 pi) would
  # miss by order one on an nll of order 300, so 1e-2 discriminates
  expect_lt(abs(-as.numeric(logLik(fit)) - cr), 1e-2)
})

test_that("anova() applies the REML design rule to non-gaussian fits", {
  dd <- rng_data(seed = 73)
  r0 <- frm(bf(yp ~ x + (1 | g)), family = poisson(), data = dd,
            REML = TRUE)
  r1 <- frm(bf(yp ~ 1 + (1 | g)), family = poisson(), data = dd,
            REML = TRUE)
  # different fixed-effect designs: the integrated likelihoods are not
  # comparable, exactly as in the gaussian case
  expect_error(anova(r0, r1), "REML")
  m0 <- frm(bf(yp ~ x + (1 | g)), family = poisson(), data = dd)
  expect_error(anova(r0, m0), "REML")
})

test_that("distributional gaussian REML is exact classical REML (gls)", {
  skip_if_not_installed("nlme")
  # sigma with its own predictor puts its coefficients into the
  # variance-parameter set; the error-contrast derivation never cares
  # how theta enters Sigma, so REML here is Patterson-Thompson
  # exactly, and gls with a varFunc under method = "REML" is the
  # authority. Only the mean coefficients are integrated out: betad
  # stays maximized, matching nlme and the double-GLM literature
  set.seed(74)
  n <- 200L
  dd <- data.frame(x = stats::rnorm(n),
                   f = factor(rep(c("a", "b"), n / 2)),
                   v = stats::runif(n, 0.5, 2))
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x,
                       ifelse(dd$f == "b", 1.6, 0.8) * exp(0.3 * dd$v))

  fi <- frm(bf(y ~ x, sigma ~ 0 + f), family = gaussian(), data = dd,
            REML = TRUE)
  gi <- nlme::gls(y ~ x, data = dd,
                  weights = nlme::varIdent(form = ~ 1 | f),
                  method = "REML")
  expect_lt(abs(as.numeric(logLik(fi)) - as.numeric(logLik(gi))), 1e-5)
  expect_lt(max(abs(unlist(fixef(fi)$mu) - stats::coef(gi))), 1e-4)

  fe <- frm(bf(y ~ x, sigma ~ v), family = gaussian(), data = dd,
            REML = TRUE)
  ge <- nlme::gls(y ~ x, data = dd,
                  weights = nlme::varExp(form = ~ v),
                  method = "REML")
  expect_lt(abs(as.numeric(logLik(fe)) - as.numeric(logLik(ge))), 1e-5)
})

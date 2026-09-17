# The scale contract, asserted THROUGH THE GENERICS.
#
# Nothing in a fitted model's output says which scale a number is on,
# and two of them silently disagree with brms: `predict()` returns the
# linear predictor where brms returns the predictive mean, and
# `summary()` prints `sigma` on its log link where brms prints the
# back-transformed value.  `?frmtmb-scales` writes the contract down;
# this file holds it.
#
# Every assertion here goes through the exported generic rather than
# through an internal, because the contract is about what a USER gets.
# Each one is paired with its inverse, the assertion that the WRONG
# scale would satisfy, so a test that stops discriminating fails
# instead of passing quietly.
#
# Tolerances.  Where a relation is an IDENTITY (the same expression on
# both sides, reached by two routes) the yardstick is machine epsilon
# in ulps.  Where it is a MEASUREMENT the yardstick is a standard error
# the run itself reports.  No absolute number appears as a tolerance.

scale_fits <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(2026)
    n <- 400L
    dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:20, 20)))
    eta <- 8 + 0.4 * dd$x + stats::rnorm(20, 0, 0.3)[dd$g]
    dd$y <- exp(stats::rnorm(n, eta, 0.4))
    dd$logy <- log(dd$y)
    cache <<- list(
      dd = dd,
      ln = frm(bf(y ~ x + (1 | g)) + lognormal(), data = dd),
      gs = frm(bf(logy ~ x + (1 | g)) + gaussian(), data = dd)
    )
    cache
  }
})

test_that("predict() reports the LINK scale and nothing else", {
  f <- scale_fits()$ln
  expect_identical(predict(f), predict(f, type = "link"))
  # the inverse: the response scale is a different number entirely, by
  # a factor the fit itself supplies
  expect_false(isTRUE(all.equal(predict(f), fitted(f))))
  ratio <- stats::median(fitted(f) / predict(f))
  expect_gt(ratio, 100)
})

test_that("predict(type = 'response') and fitted() are the same call", {
  f <- scale_fits()$ln
  # an identity: fitted() delegates to the same family mean, so the
  # yardstick is machine epsilon, not a fitted tolerance
  a <- predict(f, type = "response")
  b <- fitted(f)
  expect_lte(max(abs(a - b)), 4 * .Machine$double.eps * max(abs(b)))
})

test_that("fitted() on a lognormal is the MEAN, not the median", {
  f <- scale_fits()$ln
  mu <- predict(f, type = "link")
  sg <- sigma(f)
  # identity check first, at full precision
  expect_lte(max(abs(fitted(f) - exp(mu + sg^2 / 2))),
             8 * .Machine$double.eps * max(fitted(f)))
  # then the discrimination. The prediction's own standard error
  # cannot do it: measured on this design the mean and the median
  # sit 0.9134 of one apart. The quantity that separates them is
  # the RATIO, exp(sigma^2 / 2), whose standard error follows by
  # the delta method from the standard error of log(sigma) that
  # summary() reports.
  #
  # Read what this statistic IS before trusting it. Substituting
  # the identity, `ratio - 1` is about sigma^2/2 and `se_ratio`
  # about sigma^2 * se(log sigma), so the whole thing reduces to
  # 1 / (2 * se(log sigma)), which is about sqrt(n/2) and does NOT
  # depend on sigma. So `> 5` is, to first order, an assertion
  # that n is over 50: it is a sanity bound, not a measurement of
  # the size of the discrepancy. It still discriminates, which is
  # the job: under mutant M1 the ratio is exactly 1 and the
  # statistic is 0.
  ratio <- stats::median(fitted(f) / exp(mu))
  se_ls <- summary(f)$coefficients$sigma[1, 2]
  se_ratio <- exp(sg^2 / 2) * sg^2 * se_ls
  expect_gt((ratio - 1) / se_ratio, 5)
  # the inverse: the median does NOT satisfy the identity
  expect_false(isTRUE(all.equal(fitted(f), exp(mu))))
})

test_that("residuals() default is on the RESPONSE scale", {
  fits <- scale_fits()
  f <- fits$ln
  y <- fits$dd$y
  expect_lte(max(abs(residuals(f) + fitted(f) - y)),
             8 * .Machine$double.eps * max(abs(y)))
  # the inverse: the residual is not formed on the link scale
  expect_false(isTRUE(all.equal(residuals(f),
                                log(y) - predict(f, type = "link"))))
})

test_that("residuals(type = 'pearson') is unitless", {
  fits <- scale_fits()
  f <- fits$ln
  r <- residuals(f, type = "pearson")
  n <- length(r)
  # the sampling error of a standard deviation over n draws is about
  # 1/sqrt(2n); six of those is the bound, and it comes from n
  expect_lt(abs(stats::sd(r) - 1), 6 / sqrt(2 * n))
  # the inverse: the response-scale residual is nowhere near unit SD
  expect_gt(stats::sd(residuals(f)), 100)
})

test_that("summary() and fixef() report every dpar on its own LINK", {
  f <- scale_fits()$ln
  cf <- fixef(f)
  s <- summary(f)
  # summary() does not back-transform: it prints what fixef() holds
  expect_equal(s$coefficients$sigma[1, 1], cf$sigma[["(Intercept)"]])
  # and sigma() is the one accessor that does back-transform
  expect_lte(abs(exp(cf$sigma[["(Intercept)"]]) - sigma(f)),
             8 * .Machine$double.eps * sigma(f))
  # the inverse: the printed estimate is NOT the response-scale sigma
  expect_false(isTRUE(all.equal(s$coefficients$sigma[1, 1], sigma(f))))
  expect_lt(s$coefficients$sigma[1, 1], 0)
})

test_that("the LINK scale is the log scale, proved by the twin fit", {
  # lognormal(y) is gaussian(log y) up to a Jacobian that carries no
  # parameters, so the two fits must agree on every LINK-scale
  # quantity. The yardstick is each fit's own standard error, which is
  # what makes this a measurement rather than a chosen number.
  fits <- scale_fits()
  ln <- fits$ln
  gs <- fits$gs
  se_int <- summary(gs)$coefficients$mu[1, 2]
  expect_lt(abs(fixef(ln)$mu[["(Intercept)"]] -
                  fixef(gs)$mu[["(Intercept)"]]) / se_int, 0.01)
  # the linear predictor itself, row by row
  a <- predict(ln, type = "link")
  b <- predict(gs, type = "link")
  s <- predict(gs, type = "link", se.fit = TRUE)$se.fit
  expect_lt(max(abs(a - b)) / stats::median(s), 0.01)
  # and the random-effect covariance, which VarCorr() reports on that
  # same link scale
  va <- varcorr_matrices(ln)[[1]][1, 1]
  vb <- varcorr_matrices(gs)[[1]][1, 1]
  expect_lt(abs(va - vb) / vb, 0.01)
  # the inverse: VarCorr() is not on the response scale, where the
  # spread of the outcome is orders of magnitude larger
  expect_gt(stats::sd(fits$dd$y) / sqrt(va), 100)
})

test_that("simulate() draws on the RESPONSE scale", {
  fits <- scale_fits()
  f <- fits$ln
  y <- fits$dd$y
  sim <- simulate(f, nsim = 1, seed = 7)[[1]]
  expect_true(all(sim > 0))
  # the log-scale mean must sit within a few standard errors of the
  # observed one, and the standard error comes from the data
  se <- stats::sd(log(y)) / sqrt(length(y))
  expect_lt(abs(mean(log(sim)) - mean(log(y))) / se, 6)
  # the inverse: draws on the link scale would be near 8, not near
  # thousands
  expect_gt(stats::median(sim), 100)
})

test_that("a log-link count fit puts predict() and fitted() an exp apart", {
  set.seed(11)
  dd <- data.frame(x = stats::rnorm(300))
  dd$cnt <- stats::rpois(300, exp(0.5 + 0.4 * dd$x))
  f <- frm(bf(cnt ~ x) + poisson(), data = dd)
  expect_lte(max(abs(fitted(f) - exp(predict(f)))),
             8 * .Machine$double.eps * max(fitted(f)))
  expect_false(isTRUE(all.equal(fitted(f), predict(f))))
  # a family with no dispersion parameter reports sigma() as 1 exactly
  expect_identical(sigma(f), 1)
})

test_that("conditional_effects() is on the RESPONSE scale", {
  # The row that matters most and was missing from the page: it is
  # the one method whose scale is the OPPOSITE of predict()'s
  # default, so a porter who plots both sees three orders of
  # magnitude between them and no document explaining it.
  f <- scale_fits()$ln
  ce <- conditional_effects(f)[[1]]
  est <- ce$estimate__
  expect_true(all(est > 0))
  # it lives inside the range fitted() covers, not near the linear
  # predictor, and the yardstick is the fit's own fitted values
  fv <- fitted(f)
  expect_gt(min(est), min(fv) / 2)
  expect_lt(max(est), max(fv) * 2)
  # the inverse: it is nowhere near the link scale
  expect_gt(stats::median(est) / stats::median(predict(f)), 100)
})

test_that("the lognormal identity holds only for a constant sigma", {
  # The page used to state exp(mu + sigma^2/2) flatly. With a
  # distributional sigma there is no scalar to put in it, and
  # sigma() says so by returning NA rather than by inventing one.
  fits <- scale_fits()
  dd <- fits$dd
  f2 <- frm(bf(y ~ x + (1 | g), sigma ~ x) + lognormal(), data = dd)
  expect_warning(s <- sigma(f2), "varies by observation")
  expect_true(is.na(s))
  mu <- predict(f2, type = "link")
  sv <- predict(f2, dpar = "sigma", type = "response")
  # the per-row sigma reproduces fitted() exactly; the scalar cannot
  expect_identical(fitted(f2), exp(mu + sv^2 / 2))
  expect_true(all(is.na(exp(mu + s^2 / 2))))
})

test_that("the lognormal identity fails under truncation", {
  # fitted() is the TRUNCATED mean there, so the naive formula is not
  # merely imprecise, it is answering a different question, and how
  # far off it is is itself an identity: per row,
  # (fitted - naive) / fitted = 1 - pnorm(-a) / pnorm(sigma - a) with
  # a = (log(lb) - mu) / sigma. So the assertion is that identity, to
  # machine precision, rather than a chosen threshold on its size,
  # which depends on where the bound sits and has no single value.
  fits <- scale_fits()
  dd <- fits$dd
  # the bound is a literal: a formula is evaluated in the data, so
  # trunc(lb = lb) looks for a column called lb and fails
  lb <- 2000
  d3 <- dd[dd$y > lb, ]
  f3 <- frm(bf(y | trunc(lb = 2000) ~ x + (1 | g)) + lognormal(),
            data = d3)
  mu <- predict(f3, type = "link")
  s <- sigma(f3)
  fv <- fitted(f3)
  naive <- exp(mu + s^2 / 2)
  a <- (log(lb) - mu) / s
  closed <- 1 - stats::pnorm(-a) / stats::pnorm(s - a)
  expect_lte(max(abs((fv - naive) / fv - closed)),
             8 * .Machine$double.eps)
  # the inverse: the naive formula is not the truncated mean, by far
  # more than the precision the identity above holds to
  expect_gt(max(abs(closed)), 1e6 * .Machine$double.eps)
  expect_true(all(fv > lb))
  # and the untruncated fit on the same design still holds it exactly
  ln <- fits$ln
  expect_identical(fitted(ln),
                   exp(predict(ln, type = "link") + sigma(ln)^2 / 2))
})

# mixture() used to start every component at the quantile of the raw
# response. For a bounded mean that is a count outside the logit's
# range, so every component fell back to the link origin, where they
# are the same distribution, and the fit stayed there in silence.

test_that("a bounded-mean mixture starts from the proportion, clamped", {
  # the malingering data of Lee and Wagenmakers: eight of twenty-two
  # respondents score 45 of 45, so the two-thirds quantile of the
  # proportion is exactly 1 and only a clamp keeps the logit finite
  d <- data.frame(k = c(45, 45, 44, 45, 44, 45, 45, 45, 45, 45, 30,
                        20, 6, 44, 44, 27, 25, 17, 14, 27, 35, 30),
                  n = 45L)
  fit <- suppressWarnings(
    frm(bf(k | trials(n) ~ 1), family = mixture(beta_binomial, beta_binomial),
        data = d))
  mu <- sort(plogis(c(unname(fixef(fit)$mu1), unname(fixef(fit)$mu2))))
  expect_lt(mu[1], 0.7)
  expect_gt(mu[2], 0.9)
  # the degenerate start reached -64.26; the separated optimum is -58.71
  expect_gt(as.numeric(logLik(fit)), -59)
})

test_that("an unbounded-mean mixture keeps its quantile start", {
  set.seed(7)
  d <- data.frame(y = c(rep(0L, 30), rpois(30, 8)))
  fit <- suppressWarnings(
    frm(bf(y ~ 1), family = mixture(poisson, poisson), data = d))
  mu <- sort(exp(c(unname(fixef(fit)$mu1), unname(fixef(fit)$mu2))))
  expect_lt(mu[1], 1)
  expect_gt(mu[2], 6)
})

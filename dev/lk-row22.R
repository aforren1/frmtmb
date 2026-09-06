
test_that("row 22: the brms link roster on (0, 1) responses", {
  skip_unless_brms_fit()

  # A and B only. What these exercise is the link, not the family: each
  # row is the plainest possible model through one of the link names
  # brms accepts, so a divergence can only be the inverse link or its
  # log-odds.
  set.seed(17)
  n <- 250
  d <- data.frame(x = rnorm(n))
  d$y <- rbinom(n, 1, pnorm(0.2 + 0.8 * d$x))

  brms_lp_check(brms::bf(y ~ x), brms::bernoulli("probit"), d,
                frm(bf(y ~ x) + bernoulli("probit"), data = d))

  # brms's Stan program for probit_approx is Phi_approx(), the logistic
  # of 0.07056 x^3 + 1.5976 x, while brms's own R-side inv_link()
  # answers pnorm() for the same name. frmtmb reproduces the Stan form,
  # and this row is what says so: it would fail against the pnorm one.
  brms_lp_check(brms::bf(y ~ x), brms::bernoulli("probit_approx"), d,
                frm(bf(y ~ x) + bernoulli("probit_approx"), data = d))

  db <- data.frame(x = rnorm(n), nt = sample(3:12, n, TRUE))
  db$y <- rbinom(n, db$nt, pcauchy(0.3 + 0.7 * db$x))
  brms_lp_check(brms::bf(y | trials(nt) ~ x),
                brms::brmsfamily("binomial", link = "cauchit"), db,
                frm(bf(y | trials(nt) ~ x) +
                      stats::binomial(link = "cauchit"), data = db))

  dbe <- data.frame(x = rnorm(n))
  mu <- pnorm(0.2 + 0.5 * dbe$x)
  dbe$y <- rbeta(n, mu * 8, (1 - mu) * 8)
  brms_lp_check(brms::bf(y ~ x), brms::Beta("probit"), dbe,
                frm(bf(y ~ x) + Beta("probit"), data = dbe))
})

test_that("row 22: the brms link roster on positive-mean responses", {
  skip_unless_brms_fit()

  # Each response is drawn THROUGH its link, so the fit sits where the
  # link is well conditioned rather than where the identity would be
  # tested against a near-degenerate mean.
  set.seed(23)
  n <- 250

  dp <- data.frame(x = rnorm(n))
  dp$y <- rpois(n, (2 + 0.4 * pmax(pmin(dp$x, 3), -3))^2)
  brms_lp_check(brms::bf(y ~ x), stats::poisson(link = "sqrt"), dp,
                frm(bf(y ~ x) + stats::poisson(link = "sqrt"), data = dp))

  dn <- data.frame(x = rnorm(n))
  dn$y <- rnbinom(n, size = 3, mu = log1p(exp(1.5 + 0.5 * dn$x)))
  brms_lp_check(brms::bf(y ~ x), brms::negbinomial("softplus"), dn,
                frm(bf(y ~ x) + negbinomial("softplus"), data = dn))
  brms_lp_check(brms::bf(y ~ x), brms::negbinomial("squareplus"), dn,
                frm(bf(y ~ x) + negbinomial("squareplus"), data = dn))

  # `1/mu^2` is stats::inverse.gaussian()'s default link and brms
  # spells it the same way
  di <- data.frame(x = rnorm(n))
  di$eta <- 2 + 0.3 * pmax(pmin(di$x, 3), -3)
  di$y <- rgamma(n, 8, 8 * sqrt(di$eta))
  brms_lp_check(brms::bf(y ~ x), stats::inverse.gaussian(), di,
                frm(bf(y ~ x) + stats::inverse.gaussian(), data = di))
})

test_that("row 22c: brms 2.23.0 cannot compile its own softit program", {
  skip_unless_brms()

  # Not a divergence, and not a skip. brms emits BOTH softit helpers
  # with a vector / vector division, which Stan rejects:
  #
  #   vector softit(vector p) { return log(expm1(-p / (p - 1))); }
  #   vector inv_softit(vector y) {
  #     return log1p_exp(y) / (1 + log1p_exp(y));
  #   }
  #
  # Stan needs `./` for element-wise division of two vectors, so NO
  # brms model with this link compiles at all in 2.23.0 and there is no
  # program to check an identity against. frmtmb's softit is asserted
  # against brms's R-side inv_link() instead, in
  # test-numerical-robustness.R, where it agrees to 1.1e-16.
  #
  # This assertion is the row: when brms fixes the emitted code it
  # fails, and the identity check below it can be written for real.
  set.seed(1)
  d <- data.frame(y = rbinom(30, 1, 0.5), x = rnorm(30))
  code <- as.character(brms::stancode(
    brms::bf(y ~ x), data = d, family = brms::bernoulli(link = "softit")))
  expect_true(grepl("vector softit(vector p)", code, fixed = TRUE))
  expect_true(grepl("log(expm1(-p / (p - 1)))", code, fixed = TRUE))
  expect_true(grepl("log1p_exp(y) / (1 + log1p_exp(y))", code,
                    fixed = TRUE))
  # and neither is the element-wise form Stan would accept
  expect_false(grepl("./", code, fixed = TRUE))
})

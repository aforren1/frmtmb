
# --- conditional_effects() bands through every new link ---------------
#
# The band is built on the scale it is symmetric on and then moved to
# the response scale. Two things go wrong there and both are asserted
# here rather than assumed: an inverse link that DECREASES sends the
# lower linear predictor to the upper response bound, and a link whose
# domain does not reach an end has no bound at all on that side.
# `1/mu^2` is where this lane first met both, because
# stats::inverse.gaussian()'s default link only became fittable when the
# registry gained the name.

# one family per new link that accepts it, and the data drawn through
# that link so the fit sits where the link is well conditioned
band_case <- function(lk, seed = 5, n = 200) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  e <- 0.3 + 0.5 * d$x
  fam <- switch(lk,
    probit = , probit_approx = , cauchit = , softit = {
      lo <- frmtmb:::frmtmb_links[[lk]]
      d$y <- stats::rbinom(n, 1, lo$linkinv(e))
      bernoulli(lk)
    },
    softplus = , squareplus = , sqrt = , log1p = {
      lo <- frmtmb:::frmtmb_links[[lk]]
      # a strictly positive mean: sqrt and log1p need the predictor kept
      # off their turning points
      ep <- 1.5 + 0.4 * pmax(pmin(d$x, 3), -3)
      d$y <- stats::rpois(n, lo$linkinv(ep))
      negbinomial(lk)
    },
    `1/mu^2` = {
      ei <- 2 + 0.3 * pmax(pmin(d$x, 3), -3)
      d$y <- stats::rgamma(n, 8, 8 * sqrt(ei))
      stats::inverse.gaussian()
    },
    inverse = {
      ei <- 2 + 0.3 * pmax(pmin(d$x, 3), -3)
      d$y <- stats::rgamma(n, 8, 8 * ei)
      stats::Gamma(link = "inverse")
    },
    stop("no case for ", lk))
  list(fit = suppressWarnings(frm(bf(y ~ x), family = fam, data = d)),
       data = d)
}

test_that("conditional_effects bands are ordered through every new link", {
  for (lk in c("probit", "probit_approx", "cauchit", "softit",
               "softplus", "squareplus", "sqrt", "log1p", "1/mu^2")) {
    cs <- band_case(lk)
    expect_identical(family(cs$fit)$links$mu$name, lk,
                     label = paste(lk, "link reached the fit"))
    r <- suppressWarnings(conditional_effects(cs$fit)[[1]])
    ok <- is.finite(r$lower__) & is.finite(r$upper__)
    # NaN is never a bound: a bound the link cannot reach is NA
    expect_false(any(is.nan(r$lower__) | is.nan(r$upper__)),
                 label = paste(lk, "band has no NaN"))
    expect_true(all(r$lower__[ok] <= r$upper__[ok]),
                label = paste(lk, "band is ordered"))
    # and the estimate sits inside its own band
    expect_true(all(r$estimate__[ok] >= r$lower__[ok] &
                      r$estimate__[ok] <= r$upper__[ok]),
                label = paste(lk, "estimate inside the band"))
    expect_true(any(ok), label = paste(lk, "band is not all NA"))
  }
})

test_that("a decreasing link does not invert the band", {
  # Gamma(inverse) inverted on 100 of 100 rows before ce_band_ends():
  # 1/eta decreases, so linkinv sent the lower linear predictor to the
  # upper response bound and the two were never swapped back. This is
  # the pre-existing half of the defect - it reproduces with the link
  # registry untouched - and `1/mu^2` is the half that was new.
  for (lk in c("inverse", "1/mu^2")) {
    cs <- band_case(lk)
    r <- suppressWarnings(conditional_effects(cs$fit)[[1]])
    ok <- is.finite(r$lower__) & is.finite(r$upper__)
    expect_gt(sum(ok), 0)
    expect_equal(sum(r$lower__[ok] > r$upper__[ok]), 0,
                 label = paste(lk, "inverted rows"))
    expect_false(any(is.nan(r$lower__) | is.nan(r$upper__)))
  }
})

test_that("a bound outside the link's domain is NA, once, and named", {
  # 1/sqrt(eta) has no value at eta <= 0, so a band wide enough to reach
  # zero has no UPPER response bound - the mean runs away as eta falls
  # to zero. The bound must be NA on that side, not NaN, and not filled
  # in from the other end.
  set.seed(11)
  n <- 400
  d <- data.frame(x = stats::rnorm(n))
  mig <- 1 / sqrt(pmax(0.5 + 0.2 * d$x, 0.05))
  d$y <- stats::rgamma(n, shape = 20, rate = 20 / mig)
  fit <- suppressWarnings(frm(bf(y ~ x),
                              family = stats::inverse.gaussian(),
                              data = d))
  expect_warning(r <- conditional_effects(fit)[[1]],
                 "1/mu\\^2.*does not reach the band")
  ok <- is.finite(r$lower__) & is.finite(r$upper__)
  expect_equal(sum(r$lower__[ok] > r$upper__[ok]), 0)
  expect_false(any(is.nan(r$lower__) | is.nan(r$upper__)))
  # the missing side is the UPPER one, and the lower bound survives
  gone <- !is.finite(r$upper__)
  expect_gt(sum(gone), 0)
  expect_true(all(is.finite(r$lower__[gone])))
  expect_true(all(is.na(r$upper__[gone])))
})

test_that("ce_band_ends orders, refuses a pole and keeps the sides", {
  # the helper on its own, where the cases can be built exactly
  inc <- frmtmb:::frmtmb_links$log
  dec <- frmtmb:::frmtmb_links[["1/mu^2"]]
  # increasing: ends arrive in order and stay
  b <- frmtmb:::ce_band_ends(inc$linkinv, c(-1, 0), c(1, 2),
                             inc$mu_eta(c(0, 1)))
  expect_equal(b$lower, exp(c(-1, 0)))
  expect_equal(b$upper, exp(c(1, 2)))
  expect_identical(b$outside, 0L)
  # decreasing: the ends are swapped, not sorted
  b2 <- frmtmb:::ce_band_ends(dec$linkinv, c(1, 4), c(4, 9),
                              dec$mu_eta(c(2, 6)))
  expect_equal(b2$lower, 1 / sqrt(c(4, 9)))
  expect_equal(b2$upper, 1 / sqrt(c(1, 4)))
  expect_identical(b2$outside, 0L)
  # an end the link cannot map: NA on the side it belongs to
  b3 <- frmtmb:::ce_band_ends(dec$linkinv, c(-1), c(4), dec$mu_eta(2))
  expect_equal(b3$lower, 0.5)
  expect_true(is.na(b3$upper))
  expect_identical(b3$outside, 1L)
  # a pole BETWEEN two finite ends: neither bounds the response, so
  # both are refused. 1/eta from -1 to 1 passes through infinity.
  invl <- frmtmb:::frmtmb_links$inverse
  b4 <- frmtmb:::ce_band_ends(invl$linkinv, -1, 1, invl$mu_eta(-1))
  expect_true(is.na(b4$lower))
  expect_true(is.na(b4$upper))
  expect_identical(b4$outside, 1L)
})

test_that("predict(se.fit) reports a standard error, not an interval", {
  # H1 is a band defect and predict() has no band: it returns
  # |mu_eta| * se_eta, which is non-negative by construction and has no
  # endpoints to order. Asserted so the claim is measured, not assumed.
  for (lk in c("1/mu^2", "inverse")) {
    cs <- band_case(lk)
    nd <- data.frame(x = c(-1.5, -0.4, 0.6, 2))
    p <- stats::predict(cs$fit, newdata = nd, type = "response",
                        se.fit = TRUE)
    expect_named(p, c("fit", "se.fit"))
    expect_true(all(is.finite(p$se.fit)))
    expect_true(all(p$se.fit >= 0))
  }
})

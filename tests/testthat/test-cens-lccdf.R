# The lccdf slot: a right-censored row scored from log S directly rather
# than from log(1 - F). dev/spline-seam-proposal.md, Part 1c.
#
# The property under test is not "closer": it is that the old form has a
# region where the value is CONSTANT and the gradient is exactly zero,
# and the new one does not.

lccdf_off <- function(fam) {
  f <- frmtmb:::as_frmtmb_family(fam)
  f$lccdf <- NULL
  f
}

cens_ll <- function(fam, y, dp) {
  frmtmb:::row_lpdf(fam, y, y, dp, list(cens = rep(1, length(y))), NULL)
}

test_that("gaussian log S is exact where log(1 - F) is -Inf", {
  z <- c(1, 5, 8, 8.3, 20, 37, 100, 500)
  dp <- list(mu = rep(0, length(z)), sigma = rep(1, length(z)))
  fam <- frmtmb:::as_frmtmb_family(gaussian())
  new <- cens_ll(fam, z, dp)
  old <- cens_ll(lccdf_off(fam), z, dp)
  truth <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)

  expect_equal(new, truth, tolerance = 1e-12)
  # the old form agrees only in the representable region, and even there
  # only to 1e-8: log(1 - F) has already lost eight digits at z = 5
  expect_equal(old[z <= 5], truth[z <= 5], tolerance = 1e-8)
  expect_gt(abs(old[z == 5] - truth[z == 5]), 1e-11)
  expect_true(all(is.infinite(old[z >= 8.3])))
  expect_true(all(is.finite(new)))
  # and it is already wrong before it is infinite
  expect_gt(abs(old[z == 8] - truth[z == 8]), 0.06)
})

test_that("the censored term's GRADIENT is the failure, not its value", {
  fam <- frmtmb:::as_frmtmb_family(gaussian())
  famo <- lccdf_off(fam)
  grad <- function(f, z) {
    tp <- RTMB::MakeTape(function(p) {
      sum(frmtmb:::row_lpdf(f, z, z, list(mu = p[1], sigma = 1),
                            list(cens = 1), NULL))
    }, numeric(1))
    as.numeric(tp$jacobian(0))
  }
  for (z in c(10, 20, 40)) {
    truth <- exp(stats::dnorm(z, log = TRUE) -
                   stats::pnorm(z, lower.tail = FALSE, log.p = TRUE))
    expect_equal(grad(fam, z), truth, tolerance = 1e-8)
    expect_false(is.finite(grad(famo, z)))
  }
  # in the representable region the two agree, to the precision the old
  # form has left there
  expect_equal(grad(fam, 5), grad(famo, 5), tolerance = 1e-8)
})

test_that("weibull, exponential and lognormal log S are closed form", {
  q <- c(1, 2, 3, 4, 5, 8)
  sh <- 3
  sc <- 1 / exp(lgamma(1 + 1 / sh))
  fam <- frmtmb:::as_frmtmb_family(weibull())
  new <- cens_ll(fam, q, list(mu = rep(1, length(q)),
                              shape = rep(sh, length(q))))
  expect_equal(new, -(q / sc)^sh, tolerance = 1e-12)
  expect_true(any(is.infinite(cens_ll(lccdf_off(fam), q,
                                      list(mu = rep(1, length(q)),
                                           shape = rep(sh, length(q)))))))

  fam <- frmtmb:::as_frmtmb_family(exponential())
  expect_equal(cens_ll(fam, q, list(mu = rep(0.05, length(q)))),
               -q / 0.05, tolerance = 1e-12)

  fam <- frmtmb:::as_frmtmb_family(lognormal())
  lq <- log(q)
  expect_equal(cens_ll(fam, q, list(mu = rep(0, length(q)),
                                    sigma = rep(0.05, length(q)))),
               stats::pnorm(lq / 0.05, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-12)
})

test_that("the families that do NOT declare lccdf are the measured ones", {
  # poisson IS censored, and still cannot declare one: RTMB's ppois
  # does not tape lower.tail = FALSE, log.p = TRUE (it reaches
  # stats::ppois and errors on an advector). inverse.gaussian gains
  # nothing, because RTMBdist's upper tail is computed on the
  # probability scale and reaches -Inf at the same log S = -34 that
  # log(1 - F) does
  for (f in list(gaussian(), lognormal(), exponential(), weibull())) {
    expect_false(is.null(frmtmb:::as_frmtmb_family(f)[["lccdf"]]),
                 label = f$family)
  }
  for (f in list(poisson(), inverse.gaussian(link = "log"))) {
    expect_null(frmtmb:::as_frmtmb_family(f)[["lccdf"]], label = f$family)
  }
})

test_that("a family with lccdf alone takes right censoring and no more", {
  skip_on_cran()
  set.seed(6061)
  n <- 300
  dd <- data.frame(x = rnorm(n))
  dd$time <- rexp(n, exp(-0.5 + 0.7 * dd$x))
  dd$cens <- as.integer(dd$time > 1.5)
  dd$time <- pmin(dd$time, 1.5)

  only_s <- frmtmb_family(
    "exp_surv", dpars = "mu", links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) -log(dpars[["mu"]]) - y / dpars[["mu"]],
    lccdf = function(q, dpars, aterms) -q / dpars[["mu"]],
    init_dpars = list(mu = function(y, aterms) mean(y))
  )
  fit <- frm(bf(time | cens(cens) ~ x), dd, only_s)
  ref <- frm(bf(time | cens(cens) ~ x), dd, exponential())
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-6)

  # left censoring needs a CDF and says so
  dd2 <- dd
  dd2$cens[1:5] <- -1
  expect_error(frm(bf(time | cens(cens) ~ x), dd2, only_s),
               "lccdf")
  expect_error(frm(bf(time | trunc(lb = 0) ~ x), dd, only_s), "lccdf")
})

test_that("right censoring under trunc(ub) stays a windowed log difference", {
  skip_on_cran()
  set.seed(6062)
  n <- 200
  dd <- data.frame(x = rnorm(n))
  dd$time <- rexp(n, exp(-0.3 + 0.5 * dd$x))
  dd <- dd[dd$time < 3, ]
  dd$cens <- as.integer(dd$time > 1)
  dd$time <- pmin(dd$time, 1)
  fit <- frm(bf(time | cens(cens) + trunc(ub = 3) ~ x), dd, exponential())
  # against the same likelihood written out on the probability scale,
  # which is accurate in this (representable) region
  dp <- frmtmb:::eval_dpars(fit)[["time"]]
  mu <- dp$mu
  S <- function(q) exp(-q / mu)
  ll <- ifelse(dd$cens == 1, log(S(dd$time) - S(3)),
               -log(mu) - dd$time / mu) - log(1 - S(3))
  expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-8)
})

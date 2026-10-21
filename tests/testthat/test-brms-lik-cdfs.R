# Log-density identity against brms for the families that gained a CDF,
# each with all four censoring codes, through brms_lp_check()
# (helper-brms.R): brms's Stan program under flat priors, evaluated at
# frmtmb's estimates, equals logLik(fit), and its gradient there is
# zero, so the two likelihoods agree on the censored rows and share the
# maximum.
#
# frmtmb reads a discrete censoring or truncation bound as brms 2.23.0
# does (?frmtmb_family, "Censoring a discrete response"), so brms is
# given the same data frame frmtmb is.
#
# Stan compiles here, once per family, cached under FRMTMB_STAN_CACHE:
#   Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")

cdfs_lp_data <- function(y, seed, discrete, step = 0.5) {
  set.seed(seed)
  n <- length(y)
  d <- data.frame(y = y, x = stats::rnorm(n),
                  cc = sample(c(-1, 0, 0, 0, 1, 2), n, TRUE))
  d$y2 <- d$y + if (discrete) 2 else step
  d
}

# Stan's AD gradient of gamma_lcdf() and gamma_lccdf() in the shape is
# wrong at frmtmb's estimates while their values agree: -0.182 against
# 0.000203 from finite differences of the same log_prob, on the
# hurdle_gamma case below (dev/cdfs-dbg16.R). For a gamma CDF the
# gradient is therefore taken by central differences of brms's log
# density, which is brms's likelihood, not Stan's derivative of it.
cdfs_fd_grad <- function(r) {
  up <- rstan::unconstrain_pars(r$sf, r$pars)
  vapply(seq_along(up), function(j) {
    e <- replace(numeric(length(up)), j, 1e-6)
    (rstan::log_prob(r$sf, up + e, adjust_transform = FALSE) -
       rstan::log_prob(r$sf, up - e, adjust_transform = FALSE)) / 2e-6
  }, 0)
}

cdfs_lp_case <- function(fam_frm, fam_brms, d, discrete, trials = NULL,
                         tol_grad = 1e-3, fd = FALSE) {
  # an interval whose top was cut back to its bottom is not an interval
  d$cc[d$cc == 2 & d$y2 <= d$y] <- 0
  if (is.null(trials)) {
    f <- bf(y | cens(cc, y2) ~ x)
    bfb <- brms::bf(y | cens(cc, y2) ~ x)
  } else {
    d$tr <- trials
    f <- bf(y | trials(tr) + cens(cc, y2) ~ x)
    bfb <- brms::bf(y | trials(tr) + cens(cc, y2) ~ x)
  }
  fit <- frm(f, family = fam_frm, data = d)
  r <- brms_lp_check(bfb, fam_brms, d, fit,
                     tol_grad = if (fd) Inf else tol_grad)
  if (fd) expect_lt(max(abs(cdfs_fd_grad(r))), tol_grad)
  r
}

test_that("continuous families: censored rows are brms's likelihood", {
  skip_unless_brms_fit()
  set.seed(101)
  n <- 150
  x <- stats::rnorm(n)
  cases <- list(
    list("student", student(), brms::student(),
         1 + 0.4 * x + 0.8 * stats::rt(n, 5)),
    list("skew_normal", skew_normal(), brms::skew_normal(),
         brms::rskew_normal(n, 1, 1, 3)),
    list("exgaussian", exgaussian(), brms::exgaussian(),
         brms::rexgaussian(n, 2, 0.5, 1)),
    list("shifted_lognormal", shifted_lognormal(), brms::shifted_lognormal(),
         0.3 + stats::rlnorm(n, 0, 0.5)),
    list("von_mises", von_mises(), brms::von_mises(),
         brms::rvon_mises(n, 0.5, 3)),
    list("Gamma", Gamma(link = "log"), Gamma(link = "log"),
         stats::rgamma(n, 3, 3 / 2)),
    list("inverse.gaussian", inverse.gaussian(), stats::inverse.gaussian(),
         brms::rinv_gaussian(n, 2, 4))
  )
  for (k in cases) {
    d <- cdfs_lp_data(k[[4]], 7, FALSE)
    if (k[[1]] == "von_mises") d$y2 <- pmin(d$y + 0.4, pi)
    # brms's inv_gaussian_lccdf() has the sign of its second term
    # wrong (log(Phi(-a) + e Phi(-b)) for log(Phi(-a) - e Phi(-b))), so a
    # right-censored row is not brms's likelihood there; left and
    # interval rows read its lcdf, which is right (dev/cdfs-brms-ig.R)
    if (k[[1]] == "inverse.gaussian") d$cc[d$cc == 1] <- 0
    r <- cdfs_lp_case(k[[2]], k[[3]], d, FALSE,
                      fd = k[[1]] == "Gamma")
    expect_true(is.finite(r$lp), info = k[[1]])
  }
  # Beta's interval tops stay inside the unit interval
  d <- cdfs_lp_data(stats::rbeta(n, 2, 3), 8, FALSE)
  d$y2 <- pmin(d$y + 0.1, 0.999)
  cdfs_lp_case(Beta(), brms::Beta(), d, FALSE)
})

test_that("the asymmetric Laplace's censored rows are brms's likelihood", {
  skip_unless_brms_fit()
  set.seed(102)
  d <- cdfs_lp_data(brms::rasym_laplace(150, 1, 1, 0.3), 9, FALSE)
  # the gradient is not asserted: the likelihood has a kink at every
  # observation and nlminb stops short of the maximum on these data
  # whether or not any row is censored (dev/cdfs-dbg14.R); the value
  # identity holds at the stop
  allow_warnings(cdfs_lp_case(asym_laplace(), brms::asym_laplace(), d, FALSE,
                              tol_grad = Inf),
                 c("false convergence", "Large maximum absolute gradient"))
})

test_that("discrete families: censored rows are brms's likelihood", {
  skip_unless_brms_fit()
  set.seed(103)
  n <- 160
  cases <- list(
    list(negbinomial(), brms::negbinomial(),
         stats::rnbinom(n, mu = 3, size = 2), NULL),
    list(geometric(), brms::geometric(), stats::rnbinom(n, mu = 3, size = 1),
         NULL),
    list(binomial(), stats::binomial(), stats::rbinom(n, 12, 0.35), 12),
    list(beta_binomial(), brms::beta_binomial(),
         brms::rbeta_binomial(n, 12, 0.35, 5), 12),
    list(poisson(), stats::poisson(), stats::rpois(n, 3), NULL)
  )
  for (k in cases) {
    d <- cdfs_lp_data(k[[3]], 11, TRUE)
    if (!is.null(k[[4]])) d$y2 <- pmin(d$y2, k[[4]])
    cdfs_lp_case(k[[1]], k[[2]], d, TRUE, trials = k[[4]])
  }
})

test_that("zero-inflated and hurdle families: censored rows are brms's", {
  skip_unless_brms_fit()
  set.seed(104)
  n <- 200
  keep <- stats::rbinom(n, 1, 0.75)
  cases <- list(
    list(zero_inflated_poisson(), brms::zero_inflated_poisson(),
         stats::rpois(n, 3) * keep, TRUE, NULL),
    list(zero_inflated_negbinomial(), brms::zero_inflated_negbinomial(),
         stats::rnbinom(n, mu = 3, size = 2) * keep, TRUE, NULL),
    list(zero_inflated_binomial(), brms::zero_inflated_binomial(),
         stats::rbinom(n, 12, 0.35) * keep, TRUE, 12),
    list(zero_inflated_beta_binomial(), brms::zero_inflated_beta_binomial(),
         brms::rbeta_binomial(n, 12, 0.35, 5) * keep, TRUE, 12),
    list(hurdle_poisson(), brms::hurdle_poisson(),
         pmax(stats::rpois(n, 3), 1) * keep, TRUE, NULL),
    list(hurdle_negbinomial(), brms::hurdle_negbinomial(),
         pmax(stats::rnbinom(n, mu = 3, size = 2), 1) * keep, TRUE, NULL),
    list(hurdle_gamma(), brms::hurdle_gamma(),
         stats::rgamma(n, 3, 3 / 2) * keep, FALSE, NULL),
    list(hurdle_lognormal(), brms::hurdle_lognormal(),
         stats::rlnorm(n, 0, 0.5) * keep, FALSE, NULL),
    list(zero_inflated_beta(), brms::zero_inflated_beta(),
         stats::rbeta(n, 2, 3) * keep, FALSE, NULL)
  )
  for (k in cases) {
    d <- cdfs_lp_data(k[[3]], 12, k[[4]], step = 0.1)
    # a zero is observed; the point mass is not censored
    d$cc[d$y == 0] <- 0
    if (!is.null(k[[5]])) d$y2 <- pmin(d$y2, k[[5]])
    if (!k[[4]] && max(d$y) < 1) d$y2 <- pmin(d$y2, 0.999)
    cdfs_lp_case(k[[1]], k[[2]], d, k[[4]], trials = k[[5]],
                 fd = identical(k[[1]]$family, "hurdle_gamma"))
  }
})

test_that("trunc(): the window is brms's normalizer", {
  skip_unless_brms_fit()
  set.seed(105)
  n <- 400
  # a discrete lower bound is inclusive in both packages
  y <- stats::rnbinom(n, mu = 4, size = 2)
  d <- data.frame(y = y, x = stats::rnorm(n))
  d <- d[d$y >= 1 & d$y <= 12, ]
  fit <- frm(y | trunc(lb = 1, ub = 12) ~ x, family = negbinomial(), data = d)
  brms_lp_check(brms::bf(y | trunc(lb = 1, ub = 12) ~ x), brms::negbinomial(),
                d, fit)
  y <- 1 + 0.8 * stats::rt(n, 5)
  d <- data.frame(y = y, x = stats::rnorm(n))
  d <- d[d$y > -1 & d$y < 3, ]
  fit <- frm(y | trunc(lb = -1, ub = 3) ~ x, family = student(), data = d)
  brms_lp_check(brms::bf(y | trunc(lb = -1, ub = 3) ~ x), brms::student(),
                d, fit)
  y <- brms::rvon_mises(n, 0.5, 2)
  d <- data.frame(y = y, x = stats::rnorm(n))
  d <- d[d$y > -2 & d$y < 2.5, ]
  fit <- frm(y | trunc(lb = -2, ub = 2.5) ~ x, family = von_mises(), data = d)
  brms_lp_check(brms::bf(y | trunc(lb = -2, ub = 2.5) ~ x), brms::von_mises(),
                d, fit)
  y <- stats::rgamma(n, 3, 3 / 2) * stats::rbinom(n, 1, 0.8)
  d <- data.frame(y = y, x = stats::rnorm(n))
  d <- d[d$y < 5, ]
  fit <- frm(y | trunc(ub = 5) ~ x, family = hurdle_gamma(), data = d)
  r <- brms_lp_check(brms::bf(y | trunc(ub = 5) ~ x),
                     brms::hurdle_gamma(), d, fit, tol_grad = Inf)
  expect_lt(max(abs(cdfs_fd_grad(r))), 1e-3)
  # a zero-inflated count under a constant bound: every row's window
  # holds the point mass
  y <- stats::rpois(n, 3) * stats::rbinom(n, 1, 0.75)
  d <- data.frame(y = y, x = stats::rnorm(n))
  d <- d[d$y <= 6, ]
  fit <- frm(y | trunc(ub = 6) ~ x, family = zero_inflated_poisson(), data = d)
  brms_lp_check(brms::bf(y | trunc(ub = 6) ~ x),
                brms::zero_inflated_poisson(), d, fit)
  # censoring inside a discrete window, and trunc(lb = 0), which is none
  y <- stats::rpois(n, 4)
  d <- data.frame(y = y, x = stats::rnorm(n))
  d <- d[d$y >= 2, ]
  d$cc <- as.numeric(d$y > 6)
  d$y <- pmin(d$y, 6)
  fit <- frm(y | cens(cc) + trunc(lb = 2) ~ x, family = poisson(), data = d)
  brms_lp_check(brms::bf(y | cens(cc) + trunc(lb = 2) ~ x), stats::poisson(),
                d, fit)
  fit <- frm(y | trunc(lb = 0) ~ x, family = poisson(), data = d)
  brms_lp_check(brms::bf(y | trunc(lb = 0) ~ x), stats::poisson(), d, fit)
})

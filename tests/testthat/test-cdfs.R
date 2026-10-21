# CDFs for the families that had none: cens() and trunc() on student,
# skew_normal, exgaussian, shifted_lognormal, asym_laplace,
# zero_inflated_asym_laplace, von_mises, Gamma, Beta, binomial,
# beta_binomial, negbinomial, geometric, com_poisson, the
# zero-inflated and hurdle families, and the log survivor that poisson
# and inverse.gaussian lacked. Through frmtmb 0.72.0 every one of these
# refused cens() and trunc() with "need a family with a CDF".
#
# The references are independent of the code under test: stats, the
# closed forms, brms's R p-functions, and numerical integration of the
# family's own density.

fam_of <- function(f) frmtmb:::as_frmtmb_family(f)

# log F and log S of a family's slots at q, with the dpars given
slot_pair <- function(fam, q, dp, av = list()) {
  f <- fam_of(fam)
  list(lF = log(f$lcdf(q, dp, av)), lS = f$lccdf(q, dp, av))
}

# The hand-rolled censored log-likelihood, brms 2.23.0's on a count as
# on a continuous response: a right-censored row at y is P(Y > y) =
# 1 - F(y), and an interval is P(y < Y <= y2) = F(y2) - F(y).
hand_ll <- function(y, cc, y2, ldens, cdf) {
  ll <- ldens(y)
  r <- cc == 1
  l <- cc == -1
  i <- cc == 2
  ll[r] <- log1p(-cdf(y)[r])
  ll[l] <- log(cdf(y)[l])
  ll[i] <- log(cdf(y2)[i] - cdf(y)[i])
  ll
}

# `top` is the support's last count (a binomial's trials): a row there
# has nothing above it, so it is neither right-censored nor the bottom
# of an interval
cens_data <- function(y, seed, discrete, top = Inf) {
  set.seed(seed)
  n <- length(y)
  cc <- sample(c(-1, 0, 0, 0, 1, 2), n, TRUE)
  cc[cc %in% c(1, 2) & y >= top] <- 0
  data.frame(y = y, cc = cc, y2 = pmin(y + if (discrete) 2 else 0.5, top))
}

test_that("student's log CDF pair matches stats::pt() in both tails", {
  z <- c(-1e6, -1e3, -50, -10, -3, -1, -0.3, -1e-4, 0, 1e-4, 0.3, 1, 3, 10,
         50, 1e3, 1e6)
  for (nu in c(1.01, 2, 3.7, 4, 30, 1e3, 4e5)) {
    p <- frmtmb:::student_cdf_log(z, nu)
    expect_equal(p$lF, stats::pt(z, nu, log.p = TRUE), tolerance = 1e-12,
                 info = format(nu))
    expect_equal(p$lS, stats::pt(z, nu, lower.tail = FALSE, log.p = TRUE),
                 tolerance = 1e-12, info = format(nu))
  }
})

test_that("the incomplete gamma pair matches stats::pgamma() in both tails", {
  # to 1e6: below the mode the quadrature arm overflowed to NaN past
  # a = 12,600 until its exponent was capped
  for (a in c(0.01, 0.3, 3, 30, 300, 1e4, 1.3e4, 1e5, 1e6)) {
    x <- sort(c(a * c(1e-6, 0.1, 0.5, 0.75, 1, 1.5, 4, 20),
                a + 1 + sqrt(a) * c(-1, 0, 1, 2.5, 8)))
    x <- x[x > 0]
    p <- frmtmb:::log_igamma_two(x, a)
    rP <- stats::pgamma(x, a, log.p = TRUE)
    rQ <- stats::pgamma(x, a, lower.tail = FALSE, log.p = TRUE)
    ok <- is.finite(rP) & is.finite(rQ)
    expect_equal(p$lP[ok], rP[ok], tolerance = 1e-12, info = format(a))
    expect_equal(p$lQ[ok], rQ[ok], tolerance = 1e-12, info = format(a))
  }
})

test_that("the incomplete beta pair holds above x = 0.6, past stats", {
  # Rmpfr::pbetaI() put stats::pbeta() 7.9e-4 out at (0.99; 1e5, 7),
  # where the pair gives -970.16014318 to all ten digits
  # (dev/cdfs-acc.R). Here the closed form at a whole second shape, a
  # finite sum, is the referee: I_x(a, 1) = x^a and I_x(1, b) =
  # 1 - (1 - x)^b
  x <- c(1e-8, 0.01, 0.3, 0.5, 0.6, 0.8, 0.99, 1 - 1e-8)
  for (a in c(0.3, 7, 3e3)) {
    p <- frmtmb:::log_ibeta_two(x, 1 - x, a, 1)
    expect_equal(p$lF, a * log(x), tolerance = 1e-12, info = format(a))
    expect_equal(p$lS, log1p(-x^a), tolerance = 1e-10, info = format(a))
    q <- frmtmb:::log_ibeta_two(x, 1 - x, 1, a)
    expect_equal(q$lS, a * log1p(-x), tolerance = 1e-12, info = format(a))
  }
})

test_that("the skew-normal pair is exact at alpha = 1 to z = +-38", {
  # F(z; 1) = Phi(z)^2, so S = Phi(-z) (1 + Phi(z)); with xi = 0 and
  # omega = 1, brms's mean is delta sqrt(2 / pi) and its sd
  # sqrt(1 - 2 delta^2 / pi)
  z <- c(-38, -20, -5, -1, -0.1, 0, 0.1, 1, 5, 20, 38)
  d <- 1 / sqrt(2)
  p <- frmtmb:::skew_normal_cdf_log(z, d * sqrt(2 / pi),
                                    sqrt(1 - 2 * d^2 / pi), 1)
  expect_equal(p$lF, 2 * stats::pnorm(z, log.p = TRUE), tolerance = 1e-12)
  expect_equal(p$lS, stats::pnorm(-z, log.p = TRUE) + log1p(stats::pnorm(z)),
               tolerance = 1e-12)
  # the subtraction Stan and brms form, Phi(z) - 2 T(z, alpha), has no
  # digits left in that tail
  expect_lt(p$lF[1], -1400)
})

test_that("the von Mises pair matches brms::pvon_mises() and its own density", {
  skip_if_not_installed("brms")
  y <- c(-pi + 1e-3, -2.5, -1, 0, 0.7, 2, pi - 1e-3)
  for (k in c(0, 0.3, 3, 30, 300)) {
    for (mu in c(-2.9, 0, 1)) {
      p <- frmtmb:::von_mises_cdf_log(y, mu, k)
      expect_equal(exp(p$lF), brms::pvon_mises(y, mu, k), tolerance = 1e-9,
                   info = paste(k, mu))
      # the far side of a concentrated distribution, against the density
      # integrated with its largest value divided out
      i0e <- besselI(k, 0, expon.scaled = TRUE)
      f <- function(t) exp(k * (cos(t - mu) - 1)) / (2 * pi * i0e)
      ref <- log(stats::integrate(f, y[1], pi, rel.tol = 1e-12)$value)
      expect_equal(p$lS[1], ref, tolerance = 1e-9, info = paste(k, mu))
    }
  }
})

test_that("poisson and inverse.gaussian gained a log survivor that is exact", {
  # log(1 - F) is -Inf here for both
  expect_equal(slot_pair(poisson(), 700, list(mu = 5))$lS,
               stats::ppois(700, 5, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-12)
  expect_true(is.infinite(log1p(-stats::ppois(700, 5))))
  skip_if_not_installed("brms")
  q <- c(0.1, 1, 5, 40)
  ig <- slot_pair(inverse.gaussian(), q, list(mu = 1, shape = 2))
  ref <- brms::pinv_gaussian(q, 1, 2, lower.tail = FALSE)
  expect_equal(ig$lS[1:3], log(ref[1:3]), tolerance = 1e-10)
  expect_lt(ig$lS[4], -30)
  expect_true(is.finite(ig$lS[4]))
  # exp(2 shape / mu) overflows from shape / mu = 355; through 0.72.0
  # the CDF was Inf * 0 = NaN at every q there, so a censored or
  # truncated fit had a NaN likelihood. Against the density integrated
  q <- c(0.5, 0.9, 1, 1.1)
  F <- fam_of(inverse.gaussian())$lcdf(q, list(mu = 1, shape = 1000), list())
  ref <- vapply(q, function(u) {
    stats::integrate(function(t) {
      RTMBdist::dinvgauss(t, mean = 1, shape = 1000)
    }, 0, u, rel.tol = 1e-12)$value
  }, 0)
  expect_equal(F[-1], ref[-1], tolerance = 1e-8)
  expect_true(is.finite(F[1]) && F[1] > 0 && F[1] < 1e-100)
})

test_that("brms's p-functions agree with the families' slots", {
  skip_if_not_installed("brms")
  skip_if_not_installed("extraDistr")
  x <- c(0, 1, 2, 5, 9)
  chk <- function(fam, q, dp, av, ref) {
    expect_equal(exp(slot_pair(fam, q, dp, av)$lF), ref, tolerance = 1e-10,
                 info = fam_of(fam)$family)
  }
  chk(zero_inflated_poisson(), x, list(mu = 2, zi = 0.2), list(),
      brms::pzero_inflated_poisson(x, 2, 0.2))
  chk(zero_inflated_negbinomial(), x, list(mu = 2, shape = 3, zi = 0.2),
      list(), brms::pzero_inflated_negbinomial(x, 2, 3, 0.2))
  chk(zero_inflated_binomial(), x, list(mu = 0.3, zi = 0.2),
      list(trials = 10), brms::pzero_inflated_binomial(x, 10, 0.3, 0.2))
  chk(zero_inflated_beta_binomial(), x, list(mu = 0.3, phi = 4, zi = 0.2),
      list(trials = 10),
      brms::pzero_inflated_beta_binomial(x, 10, 0.3, 4, 0.2))
  chk(beta_binomial(), x, list(mu = 0.3, phi = 4), list(trials = 10),
      brms::pbeta_binomial(x, 10, 0.3, 4))
  chk(hurdle_poisson(), x, list(mu = 2, hu = 0.2), list(),
      brms::phurdle_poisson(x, 2, 0.2))
  chk(hurdle_negbinomial(), x, list(mu = 2, shape = 3, hu = 0.2), list(),
      brms::phurdle_negbinomial(x, 2, 3, 0.2))
  xc <- c(0, 0.3, 1, 2.5, 9)
  chk(hurdle_gamma(), xc, list(mu = 2, shape = 3, hu = 0.2), list(),
      brms::phurdle_gamma(xc, 3, 2 / 3, 0.2))
  chk(hurdle_lognormal(), xc, list(mu = 0.2, sigma = 0.7, hu = 0.2), list(),
      brms::phurdle_lognormal(xc, 0.2, 0.7, 0.2))
  xb <- c(0, 0.1, 0.5, 0.9)
  chk(zero_inflated_beta(), xb, list(mu = 0.4, phi = 5, zi = 0.2), list(),
      brms::pzero_inflated_beta(xb, 2, 3, 0.2))
  xr <- c(-4, -1, 0, 0.5, 3)
  chk(asym_laplace(), xr, list(mu = 0.3, sigma = 1.2, quantile = 0.3),
      list(), brms::pasym_laplace(xr, 0.3, 1.2, 0.3))
  chk(exgaussian(), xr, list(mu = 0.3, sigma = 1.2, beta = 0.8), list(),
      brms::pexgaussian(xr, 0.3, 1.2, 0.8))
  chk(student(), xr, list(mu = 0.3, sigma = 1.2, nu = 5), list(),
      brms::pstudent_t(xr, 5, 0.3, 1.2))
})

# One fit per family with all four censoring codes, against the
# hand-rolled likelihood at the fitted parameters. The references are
# stats, brms's p-functions, or the family's own density integrated.
cens_cases <- function() {
  set.seed(31)
  n <- 160
  x <- stats::rnorm(n)
  lin <- function(a, b) a + b * x
  list(
    student = list(
      fam = student(), y = lin(1, 0.5) + 0.8 * stats::rt(n, 5),
      ld = function(y, dp) stats::dt((y - dp$mu) / dp$sigma, dp$nu,
                                     log = TRUE) - log(dp$sigma),
      cdf = function(q, dp) stats::pt((q - dp$mu) / dp$sigma, dp$nu)),
    Gamma = list(
      fam = Gamma(link = "log"),
      y = stats::rgamma(n, 3, 3 / exp(lin(0.5, 0.3))),
      ld = function(y, dp) stats::dgamma(y, dp$shape, dp$shape / dp$mu,
                                         log = TRUE),
      cdf = function(q, dp) stats::pgamma(q, dp$shape, dp$shape / dp$mu)),
    Beta = list(
      fam = Beta(), y = stats::rbeta(n, 2 * exp(0.3 * x), 3),
      ld = function(y, dp) stats::dbeta(y, dp$mu * dp$phi,
                                        (1 - dp$mu) * dp$phi, log = TRUE),
      cdf = function(q, dp) stats::pbeta(q, dp$mu * dp$phi,
                                         (1 - dp$mu) * dp$phi)),
    negbinomial = list(
      fam = negbinomial(), y = stats::rnbinom(n, mu = exp(lin(1, 0.3)),
                                              size = 2),
      ld = function(y, dp) stats::dnbinom(y, mu = dp$mu, size = dp$shape,
                                          log = TRUE),
      cdf = function(q, dp) stats::pnbinom(q, mu = dp$mu, size = dp$shape),
      discrete = TRUE),
    geometric = list(
      fam = geometric(), y = stats::rnbinom(n, mu = exp(lin(1, 0.3)),
                                            size = 1),
      ld = function(y, dp) stats::dnbinom(y, mu = dp$mu, size = 1,
                                          log = TRUE),
      cdf = function(q, dp) stats::pnbinom(q, mu = dp$mu, size = 1),
      discrete = TRUE),
    binomial = list(
      fam = binomial(), y = stats::rbinom(n, 12, stats::plogis(0.3 * x)),
      trials = 12,
      ld = function(y, dp) stats::dbinom(y, 12, dp$mu, log = TRUE),
      cdf = function(q, dp) stats::pbinom(q, 12, dp$mu), discrete = TRUE),
    poisson = list(
      fam = poisson(), y = stats::rpois(n, exp(lin(1, 0.3))),
      ld = function(y, dp) stats::dpois(y, dp$mu, log = TRUE),
      cdf = function(q, dp) stats::ppois(q, dp$mu), discrete = TRUE),
    shifted_lognormal = list(
      fam = shifted_lognormal(), y = 0.3 + stats::rlnorm(n, 0.2 * x, 0.5),
      ld = function(y, dp) stats::dlnorm(y - dp$ndt, dp$mu, dp$sigma,
                                         log = TRUE),
      cdf = function(q, dp) stats::plnorm(q - dp$ndt, dp$mu, dp$sigma))
  )
}

test_that("each family's censored fit is its hand-rolled likelihood", {
  cs <- cens_cases()
  for (nm in names(cs)) {
    k <- cs[[nm]]
    disc <- isTRUE(k$discrete)
    d <- cens_data(k$y, 3, disc, k$trials %||% Inf)
    if (nm == "Beta") d$y2 <- pmin(d$y + 0.1, 0.999)
    # the covariate the response was drawn on, by its seed
    set.seed(31)
    d$x <- stats::rnorm(nrow(d))
    lhs <- if (!is.null(k$trials)) {
      d$tr <- k$trials
      "y | trials(tr) + cens(cc, y2)"
    } else "y | cens(cc, y2)"
    fit <- frm(stats::as.formula(paste(lhs, "~ x")), family = k$fam,
               data = d)
    dp <- frmtmb:::eval_dpars(fit)[[1]]
    ll <- hand_ll(d$y, d$cc, d$y2, function(y) k$ld(y, dp),
                  function(q) k$cdf(q, dp))
    expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-9,
                 info = nm)
    expect_equal(fit$opt$convergence, 0L, info = nm)
  }
})

test_that("the zero-inflated and hurdle fits are hand-rolled likelihoods", {
  skip_if_not_installed("brms")
  skip_if_not_installed("extraDistr")
  set.seed(41)
  n <- 200
  x <- stats::rnorm(n)
  yz <- stats::rpois(n, exp(1 + 0.3 * x)) * stats::rbinom(n, 1, 0.8)
  yh <- yz
  yg <- stats::rgamma(n, 2, 2 / exp(0.4 + 0.3 * x)) * stats::rbinom(n, 1, 0.8)
  ynb <- stats::rnbinom(n, mu = exp(1 + 0.3 * x), size = 2) *
    stats::rbinom(n, 1, 0.8)
  yb <- stats::rbinom(n, 10, stats::plogis(0.3 * x)) * stats::rbinom(n, 1, 0.8)
  ybeta <- stats::rbeta(n, 2, 3) * stats::rbinom(n, 1, 0.8)
  yln <- stats::rlnorm(n, 0.2 * x, 0.6) * stats::rbinom(n, 1, 0.8)
  # brms's d/p functions as the references; brms's beta-binomial and
  # its zero-inflated form need extraDistr
  cases <- list(
    zero_inflated_poisson = list(y = yz, disc = TRUE,
      ld = function(y, dp) brms::dzero_inflated_poisson(y, dp$mu, dp$zi,
                                                        log = TRUE),
      cdf = function(q, dp) brms::pzero_inflated_poisson(q, dp$mu, dp$zi)),
    zero_inflated_negbinomial = list(y = ynb, disc = TRUE,
      ld = function(y, dp) brms::dzero_inflated_negbinomial(
        y, dp$mu, dp$shape, dp$zi, log = TRUE),
      cdf = function(q, dp) brms::pzero_inflated_negbinomial(
        q, dp$mu, dp$shape, dp$zi)),
    zero_inflated_binomial = list(y = yb, disc = TRUE, trials = 10,
      ld = function(y, dp) brms::dzero_inflated_binomial(y, 10, dp$mu, dp$zi,
                                                         log = TRUE),
      cdf = function(q, dp) brms::pzero_inflated_binomial(q, 10, dp$mu,
                                                          dp$zi)),
    zero_inflated_beta_binomial = list(y = yb, disc = TRUE, trials = 10,
      ld = function(y, dp) brms::dzero_inflated_beta_binomial(
        y, 10, dp$mu, dp$phi, dp$zi, log = TRUE),
      cdf = function(q, dp) brms::pzero_inflated_beta_binomial(
        q, 10, dp$mu, dp$phi, dp$zi)),
    beta_binomial = list(y = yb, disc = TRUE, trials = 10,
      ld = function(y, dp) brms::dbeta_binomial(y, 10, dp$mu, dp$phi,
                                                log = TRUE),
      cdf = function(q, dp) brms::pbeta_binomial(q, 10, dp$mu, dp$phi)),
    hurdle_poisson = list(y = yh, disc = TRUE,
      ld = function(y, dp) brms::dhurdle_poisson(y, dp$mu, dp$hu, log = TRUE),
      cdf = function(q, dp) brms::phurdle_poisson(q, dp$mu, dp$hu)),
    hurdle_negbinomial = list(y = ynb, disc = TRUE,
      ld = function(y, dp) brms::dhurdle_negbinomial(y, dp$mu, dp$shape,
                                                     dp$hu, log = TRUE),
      cdf = function(q, dp) brms::phurdle_negbinomial(q, dp$mu, dp$shape,
                                                      dp$hu)),
    hurdle_gamma = list(y = yg, disc = FALSE,
      ld = function(y, dp) brms::dhurdle_gamma(y, dp$shape,
                                               dp$mu / dp$shape, dp$hu,
                                               log = TRUE),
      cdf = function(q, dp) brms::phurdle_gamma(q, dp$shape, dp$mu / dp$shape,
                                                dp$hu)),
    hurdle_lognormal = list(y = yln, disc = FALSE,
      ld = function(y, dp) brms::dhurdle_lognormal(y, dp$mu, dp$sigma,
                                                   dp$hu, log = TRUE),
      cdf = function(q, dp) brms::phurdle_lognormal(q, dp$mu, dp$sigma,
                                                    dp$hu)),
    zero_inflated_beta = list(y = ybeta, disc = FALSE,
      ld = function(y, dp) brms::dzero_inflated_beta(
        y, dp$mu * dp$phi, (1 - dp$mu) * dp$phi, dp$zi, log = TRUE),
      cdf = function(q, dp) brms::pzero_inflated_beta(
        q, dp$mu * dp$phi, (1 - dp$mu) * dp$phi, dp$zi))
  )
  for (nm in names(cases)) {
    k <- cases[[nm]]
    d <- cens_data(k$y, 5, k$disc, k$trials %||% Inf)
    d$x <- x
    if (nm == "zero_inflated_beta") d$y2 <- pmin(d$y + 0.1, 0.999)
    # a censored zero is still a zero; censor positive rows only
    d$cc[d$y == 0] <- 0
    if (!is.null(k$trials)) {
      d$tr <- k$trials
      d$y2 <- pmin(d$y2, k$trials)
    }
    f <- if (is.null(k$trials)) y | cens(cc, y2) ~ x else
      y | trials(tr) + cens(cc, y2) ~ x
    fit <- frm(f, family = brmsfamily(nm), data = d)
    dp <- frmtmb:::eval_dpars(fit)[[1]]
    ll <- hand_ll(d$y, d$cc, d$y2, function(y) k$ld(y, dp),
                  function(q) k$cdf(q, dp))
    expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-9,
                 info = nm)
  }
})

test_that("the other continuous fits are their hand-rolled likelihoods", {
  skip_if_not_installed("brms")
  set.seed(43)
  n <- 160
  x <- stats::rnorm(n)
  # the skew normal's CDF as the integral of its own density, which
  # shares no code with the angle integral under test
  sn_cdf <- function(q, dp) {
    m <- rep_len(dp$mu, length(q))
    s <- rep_len(dp$sigma, length(q))
    a <- rep_len(dp$alpha, length(q))
    vapply(seq_along(q), function(i) {
      stats::integrate(function(t) {
        RTMBdist::dskewnorm2(t, m[i], s[i], a[i])
      }, -Inf, q[i], rel.tol = 1e-12)$value
    }, 0)
  }
  cases <- list(
    skew_normal = list(
      y = brms::rskew_normal(n, 1 + 0.5 * x, 1, 4),
      ld = function(y, dp) brms::dskew_normal(y, dp$mu, dp$sigma, dp$alpha,
                                              log = TRUE),
      cdf = sn_cdf),
    exgaussian = list(
      y = brms::rexgaussian(n, 2 + 0.5 * x, 0.5, 1),
      ld = function(y, dp) brms::dexgaussian(y, dp$mu, dp$sigma, dp$beta,
                                             log = TRUE),
      cdf = function(q, dp) brms::pexgaussian(q, dp$mu, dp$sigma, dp$beta)),
    asym_laplace = list(
      y = brms::rasym_laplace(n, 1 + 0.5 * x, 1, 0.3),
      ld = function(y, dp) brms::dasym_laplace(y, dp$mu, dp$sigma,
                                               dp$quantile, log = TRUE),
      cdf = function(q, dp) brms::pasym_laplace(q, dp$mu, dp$sigma,
                                                dp$quantile)),
    von_mises = list(
      y = brms::rvon_mises(n, 2 * atan(0.4 * x), 3),
      ld = function(y, dp) brms::dvon_mises(y, dp$mu, dp$kappa, log = TRUE),
      cdf = function(q, dp) brms::pvon_mises(q, dp$mu, dp$kappa)),
    inverse.gaussian = list(
      y = brms::rinv_gaussian(n, exp(0.5 + 0.2 * x), 4),
      ld = function(y, dp) brms::dinv_gaussian(y, dp$mu, dp$shape, log = TRUE),
      cdf = function(q, dp) brms::pinv_gaussian(q, dp$mu, dp$shape))
  )
  for (nm in names(cases)) {
    k <- cases[[nm]]
    d <- cens_data(k$y, 6, FALSE)
    d$x <- x
    if (nm == "von_mises") d$y2 <- pmin(d$y + 0.4, pi)
    fam <- if (nm == "inverse.gaussian") {
      inverse.gaussian(link = "log")
    } else brmsfamily(nm)
    # the asymmetric Laplace likelihood has a kink at every observation,
    # and its fit stops at nlminb's false convergence on these data with
    # no censoring at all, on the 0.72.0 build as on this one
    # (dev/cdfs-dbg14.R); the identity below holds at the stop
    ok <- if (nm == "asym_laplace") {
      c("false convergence", "Large maximum absolute gradient")
    } else character()
    fit <- allow_warnings(frm(y | cens(cc, y2) ~ x, family = fam, data = d),
                          ok)
    dp <- frmtmb:::eval_dpars(fit)[[1]]
    ll <- hand_ll(d$y, d$cc, d$y2, function(y) k$ld(y, dp),
                  function(q) k$cdf(q, dp))
    expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-8,
                 info = nm)
  }
  # compois: its CDF as the sum of frmtmb's own mass function, which is
  # RTMB::dcompois2() and shares nothing with the recursion under test
  yc <- stats::rpois(n, exp(1 + 0.3 * x))
  d <- cens_data(yc, 7, TRUE)
  d$x <- x
  fit <- frm(y | cens(cc, y2) ~ x, family = brmsfamily("com_poisson"),
             data = d)
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  nu <- rep_len(dp$nu, n)
  cdf <- function(q) {
    vapply(seq_along(q), function(i) {
      if (q[i] < 0) return(0)
      sum(exp(RTMB::dcompois2(0:q[i], dp$mu[i], nu[i], log = TRUE)))
    }, 0)
  }
  ll <- hand_ll(d$y, d$cc, d$y2,
                function(y) RTMB::dcompois2(y, dp$mu, nu, log = TRUE),
                cdf)
  expect_equal(as.numeric(logLik(fit)), sum(ll), tolerance = 1e-8)
})

test_that("trunc() fits, and the truncated means are brute-force sums", {
  # brms 2.23.0 has a truncated mean for student, Gamma, binomial,
  # negbinomial and geometric; the lower bound of a count is inclusive,
  # as the likelihood reads it
  mu <- 3.2
  sh <- 1.7
  lb <- 2
  ub <- 9
  k <- lb:ub
  w <- stats::dnbinom(k, mu = mu, size = sh)
  expect_equal(frmtmb:::nbinom_trunc_mean(mu, sh, lb, ub), sum(k * w) / sum(w),
               tolerance = 1e-12)
  w <- stats::dbinom(k, 12, 0.4)
  expect_equal(frmtmb:::binomial_trunc_mean(0.4, 12, lb, ub),
               sum(k * w) / sum(w), tolerance = 1e-12)
  f <- function(t) t * stats::dt((t - 1) / 2, 5) / 2
  m <- stats::integrate(f, -1, 4)$value /
    (stats::pt((4 - 1) / 2, 5) - stats::pt((-1 - 1) / 2, 5))
  expect_equal(frmtmb:::student_trunc_mean(1, 2, 5, -1, 4), m,
               tolerance = 1e-10)
  g <- function(t) t * stats::dgamma(t, 3, 3 / 2)
  m <- stats::integrate(g, 0.5, 4)$value /
    diff(stats::pgamma(c(0.5, 4), 3, 3 / 2))
  expect_equal(frmtmb:::gamma_trunc_mean(2, 3, 0.5, 4), m, tolerance = 1e-10)

  set.seed(7)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rnbinom(300, mu = exp(1 + 0.3 * d$x), size = 2)
  d <- d[d$y >= 1 & d$y <= 12, ]
  fit <- frm(y | trunc(lb = 1, ub = 12) ~ x, family = negbinomial(), data = d)
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  hand <- stats::dnbinom(d$y, mu = dp$mu, size = dp$shape, log = TRUE) -
    log(stats::pnbinom(12, mu = dp$mu, size = dp$shape) -
          stats::pnbinom(0, mu = dp$mu, size = dp$shape))
  expect_equal(as.numeric(logLik(fit)), sum(hand), tolerance = 1e-9)
  fv <- fitted(fit)
  ref <- frmtmb:::nbinom_trunc_mean(dp$mu, dp$shape, 1, 12)
  expect_equal(unname(fv[, 1]), unname(ref), tolerance = 1e-10)
})

test_that("a family brms refuses is refused, with the registry's reason", {
  d <- data.frame(y = rep(0:1, 20), x = stats::rnorm(40), cc = 0)
  expect_error(frm(y | cens(cc) ~ x, family = bernoulli(), data = d),
               "need a family with a CDF")
  expect_equal(frm_compat("cens()", "bernoulli")$status, "refused")
  expect_match(frm_compat("trunc()", "xbeta")$note, "brms 2.23.0 refuses")
})

test_that("the censored objective's gradient is the finite difference's", {
  set.seed(9)
  d <- data.frame(x = stats::rnorm(80))
  d$y <- 1 + 0.5 * d$x + stats::rt(80, 4)
  d$cc <- sample(c(-1, 0, 1), 80, TRUE)
  fit <- frm(y | cens(cc) ~ x, family = student(), data = d)
  p <- fit$obj$par + 0.05
  g <- fit$obj$gr(p)
  fd <- vapply(seq_along(p), function(j) {
    h <- 1e-5 * max(1, abs(p[j]))
    e <- replace(numeric(length(p)), j, h)
    (fit$obj$fn(p + e) - fit$obj$fn(p - e)) / (2 * h)
  }, 0)
  expect_equal(as.numeric(g), fd, tolerance = 1e-6)
})

test_that("the CDF at z = 0 has the density as its slope", {
  # every form of the student and skew-normal tails reads z through
  # z^2 or |z|; the slope at zero comes from the forms written for it
  f <- function(p) {
    lp <- frmtmb:::student_cdf_log(0.3 - p[1], 4)
    lp$lF
  }
  tp <- RTMB::MakeTape(f, 0.3)
  expect_equal(as.numeric(tp$jacobian(0.3)),
               -stats::dt(0, 4) / stats::pt(0, 4), tolerance = 1e-10)
  g <- function(p) frmtmb:::skew_normal_cdf_log(0, p[1], 1, 2)$lS
  tg <- RTMB::MakeTape(g, 0)
  num <- (frmtmb:::skew_normal_cdf_log(0, 1e-6, 1, 2)$lS -
            frmtmb:::skew_normal_cdf_log(0, -1e-6, 1, 2)$lS) / 2e-6
  expect_equal(as.numeric(tg$jacobian(0)), num, tolerance = 1e-6)
})

test_that("a censored Gamma model with a random intercept fits", {
  # the Laplace approximation needs third derivatives, which
  # RTMB::pgamma() does not have at about half of all points
  set.seed(13)
  g <- factor(rep(1:15, each = 12))
  u <- stats::rnorm(15, 0, 0.4)[g]
  d <- data.frame(g = g, x = stats::rnorm(180))
  d$y <- stats::rgamma(180, 3, 3 / exp(0.5 + 0.3 * d$x + u))
  d$cc <- as.numeric(d$y > stats::quantile(d$y, 0.8))
  d$y[d$cc == 1] <- stats::quantile(d$y, 0.8)
  fit <- frm(y | cens(cc) ~ x + (1 | g), family = Gamma(link = "log"),
             data = d)
  expect_equal(fit$opt$convergence, 0L)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

test_that("a censored Gamma at shape 1e5 reaches the independent maximum", {
  skip_on_cran()
  # shape 1e5 is a CV of 0.3 percent; past a = 12,600 the incomplete
  # gamma below the mode was NaN and the fit stopped 20 below the top
  set.seed(3)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rgamma(200, shape = 1e5, rate = 1e5 / exp(1 + 0.1 * d$x))
  qq <- stats::quantile(d$y, 0.8)
  d$cc <- as.integer(d$y > qq)
  d$y[d$cc == 1] <- qq
  fit <- frm(y | cens(cc) ~ x, family = Gamma(link = "log"), data = d)
  nll <- function(p) {
    mu <- exp(p[1] + p[2] * d$x)
    sh <- exp(p[3])
    -sum(ifelse(d$cc == 1,
                stats::pgamma(d$y, sh, sh / mu, lower.tail = FALSE,
                              log.p = TRUE),
                stats::dgamma(d$y, sh, sh / mu, log = TRUE)))
  }
  ind <- stats::nlminb(c(1, 0.1, log(1e5)), nll)
  expect_equal(fit$opt$convergence, 0L)
  expect_equal(as.numeric(logLik(fit)), -ind$objective, tolerance = 1e-9)
})

test_that("poisson right censoring at counts near 1e5 is exact", {
  skip_on_cran()
  # the survivor is the lower incomplete gamma at a = q + 1, so a count
  # past 12,600 reached the same NaN
  expect_equal(frmtmb:::poisson_cdf_log(1e5, 7e4)$lS,
               stats::ppois(1e5, 7e4, lower.tail = FALSE, log.p = TRUE),
               tolerance = 1e-12)
  set.seed(11)
  d <- data.frame(x = stats::rnorm(300))
  d$y <- stats::rpois(300, exp(11 + 0.3 * d$x))
  d$cc <- as.numeric(d$y > 75000)
  d$y[d$cc == 1] <- 75000
  fit <- frm(y | cens(cc) ~ x, family = poisson(), data = d)
  nll <- function(p) {
    mu <- exp(p[1] + p[2] * d$x)
    -sum(ifelse(d$cc == 1,
                stats::ppois(d$y, mu, lower.tail = FALSE, log.p = TRUE),
                stats::dpois(d$y, mu, log = TRUE)))
  }
  ind <- stats::nlminb(c(11, 0.3), nll)
  expect_equal(fit$opt$convergence, 0L)
  expect_equal(as.numeric(logLik(fit)), -ind$objective, tolerance = 1e-9)
})

# The CDF of each zero-inflated and hurdle family by hand, at q (data)
# and the dpars of a fit, from stats, extraDistr and brms.
zi_hurdle_cdf <- function(nm, q, dp, tr = NULL) {
  pm <- function(p) ifelse(q < 0, 0, p)
  zi <- dp$zi
  hu <- dp$hu
  switch(nm,
    zero_inflated_poisson = pm(zi + (1 - zi) * stats::ppois(q, dp$mu)),
    zero_inflated_negbinomial =
      pm(zi + (1 - zi) * stats::pnbinom(q, size = dp$shape, mu = dp$mu)),
    zero_inflated_binomial = pm(zi + (1 - zi) * stats::pbinom(q, tr, dp$mu)),
    zero_inflated_beta_binomial =
      pm(zi + (1 - zi) * extraDistr::pbbinom(q, tr, dp$mu * dp$phi,
                                             (1 - dp$mu) * dp$phi)),
    zero_inflated_beta =
      pm(zi + (1 - zi) * stats::pbeta(q, dp$mu * dp$phi,
                                      (1 - dp$mu) * dp$phi)),
    zero_inflated_asym_laplace =
      zi * (q >= 0) + (1 - zi) * brms::pasym_laplace(q, dp$mu, dp$sigma,
                                                     dp$quantile),
    hurdle_poisson = {
      p0 <- stats::dpois(0, dp$mu)
      pm(hu + (1 - hu) * (stats::ppois(q, dp$mu) - p0) / (1 - p0))
    },
    hurdle_negbinomial = {
      p0 <- stats::dnbinom(0, size = dp$shape, mu = dp$mu)
      pm(hu + (1 - hu) * (stats::pnbinom(q, size = dp$shape, mu = dp$mu) -
                            p0) / (1 - p0))
    },
    hurdle_gamma =
      pm(hu + (1 - hu) * stats::pgamma(q, dp$shape, dp$shape / dp$mu)),
    hurdle_lognormal =
      pm(hu + (1 - hu) * stats::plnorm(q, dp$mu, dp$sigma)))
}

test_that("a constant trunc() bound gives every row the point mass", {
  skip_if_not_installed("extraDistr")
  skip_if_not_installed("brms")
  # zi_cdf_log() sized its rows by length(q), and a constant bound is one
  # q, so the point mass reached row 1 alone and the optimizer ran zi to
  # one. Each family: the bound as a constant, the same bound as a
  # column, and the bound varying by row, each objective at the same
  # parameters against the hand-built likelihood there.
  set.seed(51)
  n <- 120
  x <- stats::rnorm(n)
  zi <- stats::runif(n) < 0.25
  yc <- stats::rpois(n, exp(1 + 0.3 * x))
  ybin <- stats::rbinom(n, 10, stats::plogis(0.3 * x))
  ycont <- stats::rgamma(n, 3, 3 / exp(0.2 * x))
  cases <- list(
    zero_inflated_poisson = list(
      fam = zero_inflated_poisson(), y = ifelse(zi, 0, yc), ub = 30,
      b = c(1, 0.3), dp = list(zi = 0.25)),
    zero_inflated_negbinomial = list(
      fam = zero_inflated_negbinomial(), y = ifelse(zi, 0, yc), ub = 30,
      b = c(1, 0.3), dp = list(shape = 4, zi = 0.25)),
    zero_inflated_binomial = list(
      fam = zero_inflated_binomial(), y = ifelse(zi, 0, ybin), ub = 9,
      tr = 10, b = c(0, 0.3), dp = list(zi = 0.25)),
    zero_inflated_beta_binomial = list(
      fam = zero_inflated_beta_binomial(), y = ifelse(zi, 0, ybin), ub = 9,
      tr = 10, b = c(0, 0.3), dp = list(phi = 6, zi = 0.25)),
    zero_inflated_beta = list(
      fam = zero_inflated_beta(), y = ifelse(zi, 0, stats::rbeta(n, 2, 3)),
      ub = 0.95, b = c(-0.4, 0.2), dp = list(phi = 5, zi = 0.25)),
    zero_inflated_asym_laplace = list(
      fam = zero_inflated_asym_laplace(),
      y = ifelse(zi, 0, 1 + 0.3 * x + stats::rnorm(n)), ub = 6,
      b = c(1, 0.3), dp = list(sigma = 0.8, quantile = 0.4, zi = 0.25)),
    hurdle_poisson = list(
      fam = hurdle_poisson(), y = ifelse(zi, 0, pmax(yc, 1)), ub = 30,
      b = c(1, 0.3), dp = list(hu = 0.25)),
    hurdle_negbinomial = list(
      fam = hurdle_negbinomial(), y = ifelse(zi, 0, pmax(yc, 1)), ub = 30,
      b = c(1, 0.3), dp = list(shape = 4, hu = 0.25)),
    hurdle_gamma = list(
      fam = hurdle_gamma(), y = ifelse(zi, 0, ycont), ub = 30,
      b = c(0, 0.2), dp = list(shape = 3, hu = 0.25)),
    hurdle_lognormal = list(
      fam = hurdle_lognormal(), y = ifelse(zi, 0, ycont), ub = 30,
      b = c(0, 0.2), dp = list(sigma = 0.6, hu = 0.25)))
  for (nm in names(cases)) {
    k <- cases[[nm]]
    fam <- fam_of(k$fam)
    # the constant bound has to hold every row: brms refuses a response
    # outside the window, and so does frmtmb
    y <- pmin(k$y, k$ub)
    disc <- identical(fam$type, "discrete")
    d <- data.frame(x = x, y = y, tr = k$tr %||% 1)
    d$ubc <- k$ub
    # per row: an upper bound above each response, and on a count a
    # lower bound of 0 (absent) or 1 where the response allows it
    d$ubr <- if (disc) y + sample(1:3, n, TRUE) else y + stats::runif(n)
    d$lbr <- if (disc) ifelse(y >= 1, sample(0:1, n, TRUE), 0) else -Inf
    tri <- if (!is.null(k$tr)) "trials(tr) + " else ""
    obj <- function(term) {
      f <- stats::as.formula(paste0("y | ", tri, term, " ~ x"))
      frm(f, family = k$fam, data = d, dry_run = "objective")$obj
    }
    # the parameters through the family's own links, mu's two
    # coefficients first and the other dpars' intercepts in order
    others <- setdiff(fam$dpars, "mu")
    p <- c(k$b, vapply(others, function(o) fam$links[[o]]$linkfun(k$dp[[o]]),
                       0))
    dp <- c(list(mu = fam$links$mu$linkinv(k$b[1] + k$b[2] * x)), k$dp)
    hand <- function(lb, ub) {
      ld <- fam$lpdf(d$y, dp, list(trials = d$tr))
      lo <- if (disc) lb - 1 else lb
      Z <- zi_hurdle_cdf(nm, ub, dp, d$tr) -
        ifelse(is.finite(lo) & (!disc | lo >= 0),
               zi_hurdle_cdf(nm, lo, dp, d$tr), 0)
      sum(ld - log(Z))
    }
    ref <- hand(rep(-Inf, n), rep(k$ub, n))
    expect_equal(-obj(paste0("trunc(ub = ", k$ub, ")"))$fn(p), ref,
                 tolerance = 1e-10, info = paste(nm, "constant"))
    expect_equal(-obj("trunc(ub = ubc)")$fn(p), ref, tolerance = 1e-10,
                 info = paste(nm, "column"))
    expect_equal(-obj("trunc(lb = lbr, ub = ubr)")$fn(p),
                 hand(d$lbr, d$ubr), tolerance = 1e-10,
                 info = paste(nm, "per row"))
  }
})

test_that("the ex-gaussian density is the closed form past sigma / beta = 20", {
  skip_if_not_installed("Rmpfr")
  # RTMBdist::dexgauss() switches to a gaussian approximation near
  # u = 20, where it is 0.16 wrong in the log density
  fam <- fam_of(exgaussian())
  y <- c(-3, -1, 0, 0.5, 2, 3)
  for (u in c(5, 20, 100, 1e3)) {
    sigma <- 1
    beta <- sigma / u
    mu <- 0.2
    got <- fam$lpdf(y, list(mu = mu, sigma = sigma, beta = beta), list())
    m <- function(v) Rmpfr::mpfr(v, 300)
    z <- (m(y) - mu + m(beta)) / sigma
    uu <- m(sigma) / m(beta)
    ref <- log(exp(uu * uu / 2 - uu * z) * Rmpfr::pnorm(z - uu)) - log(m(beta))
    expect_equal(got, as.numeric(ref), tolerance = 1e-12, info = format(u))
  }
})

test_that("the ex-gaussian CDF's derivatives hold at sigma / beta of 1e3 and 1e4", {
  skip_if_not_installed("Rmpfr")
  # lx = u^2 / 2 - u z + log Phi(z - u) cancelled to O(1 / u) from terms
  # of size u^2, so the value held and the derivatives did not
  # p is 400-bit, so a step of 1e-12 is not rounded to the doubles
  ref_lS <- function(p, q) {
    mu <- p[1]
    sigma <- exp(p[2])
    beta <- exp(p[3])
    u <- sigma / beta
    z <- (Rmpfr::mpfr(q, 400) - mu + beta) / sigma
    log(Rmpfr::pnorm(-z) + exp(u * u / 2 - u * z) * Rmpfr::pnorm(z - u))
  }
  f <- function(p) frmtmb:::exgaussian_cdf_log(-1.5, p[1], exp(p[2]),
                                               exp(p[3]))$lS
  for (lb in c(log(1e-3), log(1e-4))) {
    p <- c(-2, 0, lb)
    ad <- as.numeric(RTMB::MakeTape(f, p)$jacobian(p))
    h <- Rmpfr::mpfr(1e-12, 400)
    pm <- Rmpfr::mpfr(p, 400)
    fd <- vapply(1:3, function(j) {
      e <- pm * 0
      e[j] <- h
      as.numeric((ref_lS(pm + e, -1.5) - ref_lS(pm - e, -1.5)) / (2 * h))
    }, 0)
    expect_equal(ad, fd, tolerance = 1e-8, info = format(exp(-lb)))
  }
})

test_that("compois's survivor is the upper sum where the complement is lost", {
  # the complement of a CDF near one was -690.8 with a zero gradient
  lse <- function(v) {
    mx <- max(v)
    mx + log(sum(exp(v - mx)))
  }
  for (pr in list(c(0.1, 0.3, 10), c(3, 1, 48), c(30, 3, 60))) {
    lm <- RTMB::dcompois2(0:3000, pr[1], pr[2], log = TRUE)
    ref <- lse(lm[(pr[3] + 2):3001])
    got <- frmtmb:::compois_cdf_log(pr[3], pr[1], pr[2])$lS
    expect_equal(got, ref, tolerance = 1e-12, info = paste(pr, collapse = " "))
  }
})

test_that("a hurdle's F(0) is hu to the last digit at small hu", {
  # as the complement of the survivor it was 7.8e-6 relative at 1e-12
  for (hu in c(1e-12, 1e-6, 0.3)) {
    g <- list(l = log(hu), l1m = log1p(-hu))
    p <- frmtmb:::hurdle_cdf_log(c(0, 3), frmtmb:::poisson_cdf_log(c(0, 3), 2),
                                 g, log(-expm1(-2)))
    expect_equal(p$lF[1], log(hu), tolerance = 1e-14, info = format(hu))
    p0 <- stats::dpois(0, 2)
    F3 <- hu + (1 - hu) * (stats::ppois(3, 2) - p0) / (1 - p0)
    expect_equal(p$lF[2], log(F3), tolerance = 1e-14, info = format(hu))
  }
})

test_that("a left-censored row below 1e-308 keeps its log probability", {
  # the probability-scale lcdf ends at 1e-308; cdf_log does not
  fam <- fam_of(gaussian())
  z <- c(-40, -100, -1000)
  ll <- frmtmb:::row_lpdf(fam, z, z, list(mu = rep(0, 3), sigma = rep(1, 3)),
                          list(cens = rep(-1, 3)), NULL)
  expect_equal(ll, stats::pnorm(z, log.p = TRUE), tolerance = 1e-12)
  famg <- fam_of(Gamma(link = "log"))
  q <- c(1e-4, 1e-6)
  llg <- frmtmb:::row_lpdf(famg, q, q, list(mu = rep(1, 2), shape = rep(80, 2)),
                           list(cens = rep(-1, 2)), NULL)
  expect_equal(llg, stats::pgamma(q, 80, 80, log.p = TRUE), tolerance = 1e-12)
})

test_that("a short von Mises arc next to +-pi keeps its digits", {
  # (-pi, -pi + 1e-9] was a difference of two tails that nearly agree
  # and kept 4e-8 of its relative accuracy
  for (k in c(1e-3, 0.5, 5, 50)) {
    for (mu in c(-3.1, 0, 3.1)) {
      lf <- frmtmb:::von_mises_cdf_log(c(-pi + 1e-9, pi - 1e-9), mu, k)
      dens <- function(t) exp(k * (cos(t - mu) - 1)) /
        (2 * pi * besselI(k, 0, expon.scaled = TRUE))
      # on an arc of 1e-9 the density is constant to 1e-18 relative;
      # the arcs are the lengths the doubles hold, 1e-9 to 8e-8
      L <- c((-pi + 1e-9) + pi, pi - (pi - 1e-9))
      ref <- log(L * dens(c(-pi + L[1] / 2, pi - L[2] / 2)))
      expect_equal(lf$lF[1], ref[1], tolerance = 1e-12,
                   info = paste(k, mu))
      expect_equal(lf$lS[2], ref[2], tolerance = 1e-12,
                   info = paste(k, mu))
    }
  }
})

test_that("a count's cheap forms keep both tails exact", {
  # poisson takes each row's form by its count: closed form at 0,
  # finite sums of the masses to 30, the incomplete gamma above; the
  # negative binomial the closed form at 0
  q <- c(0, 1, 3, 10, 30, 31, 60)
  for (mu in c(1e-3, 0.5, 3, 25, 200)) {
    p <- frmtmb:::poisson_cdf_log(q, mu)
    ok <- is.finite(stats::ppois(q, mu, log.p = TRUE))
    expect_equal(p$lF[ok], stats::ppois(q, mu, log.p = TRUE)[ok],
                 tolerance = 1e-12, info = format(mu))
    expect_equal(p$lS, stats::ppois(q, mu, lower.tail = FALSE, log.p = TRUE),
                 tolerance = 1e-12, info = format(mu))
  }
  for (pr in list(c(1e8, 1e-8), c(1e6, 1e-8), c(0.5, 3), c(2, 1e4))) {
    p <- frmtmb:::nbinom_cdf_log(0, log(pr[2]), log(pr[1]))
    expect_equal(p$lF, stats::dnbinom(0, size = pr[1], mu = pr[2], log = TRUE),
                 tolerance = 1e-12, info = paste(pr, collapse = " "))
    expect_equal(p$lS, stats::pnbinom(0, size = pr[1], mu = pr[2],
                                      lower.tail = FALSE, log.p = TRUE),
                 tolerance = 1e-12, info = paste(pr, collapse = " "))
  }
})

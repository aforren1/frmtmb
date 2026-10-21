# Every censored row and every trunc() window is the mass of a window
# (lo, hi], formed on the log scale from both tails of the family's CDF
# (row_lpdf(), log_window()). Through the first build of the families'
# CDFs core still formed an interval as log(F2 - F1) and the trunc(lb = )
# normalizer as log(1 - F(lb)) on the probability scale: both rounded to
# nothing in an upper tail, and at a hurdle's point mass.
#
# The references are independent: Rmpfr quadrature of the density,
# stats' own log tails, and finite differences of the objective.

win_fam <- function(f) frmtmb:::as_frmtmb_family(f)

# the objective's parameter vector at given dpars, through the family's
# own links, for an intercept-only bf() with every dpar ~ 1
win_par <- function(fam, dp) {
  f <- win_fam(fam)
  vapply(f$dpars, function(nm) f$links[[nm]]$linkfun(dp[[nm]]), 0)
}

win_obj <- function(fam, d, lhs) {
  f <- win_fam(fam)
  dpn <- setdiff(f$dpars, "mu")
  fo <- do.call(bf, c(list(stats::as.formula(paste(lhs, "~ 1"))),
                      lapply(dpn, function(dp) {
                        stats::as.formula(paste(dp, "~ 1"))
                      })))
  frm(fo, family = fam, data = d, dry_run = "objective")$obj
}

# the row's log mass by 300-bit Romberg quadrature of a density written
# in Rmpfr
mpfr_mass <- function(dens, a, b) {
  r <- Rmpfr::integrateR(dens, Rmpfr::mpfr(a, 300), Rmpfr::mpfr(b, 300),
                         rel.tol = 1e-25, verbose = FALSE)
  as.numeric(log(r$value))
}

test_that("an interval in the upper tail keeps its digits", {
  skip_if_not_installed("Rmpfr")
  # skew_normal at alpha = -10: F(y) rounds to one from y = 3, where
  # log(F2 - F1) was 0.21 off
  dp <- list(mu = 0, sigma = 1, alpha = -10)
  delta <- dp$alpha / sqrt(1 + dp$alpha^2)
  omega <- dp$sigma / sqrt(1 - 2 * delta^2 / pi)
  xi <- dp$mu - omega * delta * sqrt(2 / pi)
  dsn_m <- function(t) {
    z <- (t - xi) / omega
    2 / omega * Rmpfr::dnorm(z) * Rmpfr::pnorm(dp$alpha * z)
  }
  y <- c(1, 2, 2.6, 3.0, 3.3)
  y2 <- y + 0.7
  got <- frmtmb:::row_lpdf(win_fam(skew_normal()), y, y,
                           lapply(dp, rep, length(y)),
                           list(cens = rep(2, length(y)), cens_y2 = y2), NULL)
  ref <- mapply(function(a, b) mpfr_mass(dsn_m, a, b), y, y2)
  expect_equal(got, ref, tolerance = 1e-11)
  # gaussian, whose two log tails are stats' own
  g <- frmtmb:::row_lpdf(win_fam(gaussian()), c(9, 30), c(9, 30),
                         list(mu = c(0, 0), sigma = c(1, 1)),
                         list(cens = c(2, 2), cens_y2 = c(9.5, 31)), NULL)
  s1 <- stats::pnorm(c(9, 30), lower.tail = FALSE, log.p = TRUE)
  s2 <- stats::pnorm(c(9.5, 31), lower.tail = FALSE, log.p = TRUE)
  expect_equal(g, s1 + log(-expm1(s2 - s1)), tolerance = 1e-13)
})

test_that("a von Mises interval on the far side keeps its digits", {
  skip_if_not_installed("Rmpfr")
  k <- 50
  mu <- -2
  i0 <- Rmpfr::integrateR(function(t) exp(k * (cos(t - mu) - 1)),
                          Rmpfr::mpfr(-pi, 300), Rmpfr::mpfr(pi, 300),
                          rel.tol = 1e-25)$value
  dvm <- function(t) exp(k * (cos(t - mu) - 1)) / i0
  y <- c(0.6, 1.0, 1.154)
  y2 <- y + 0.4
  got <- frmtmb:::row_lpdf(win_fam(von_mises()), y, y,
                           list(mu = rep(mu, 3), kappa = rep(k, 3)),
                           list(cens = rep(2, 3), cens_y2 = y2), NULL)
  ref <- mapply(function(a, b) mpfr_mass(dvm, a, b), y, y2)
  expect_equal(got, ref, tolerance = 1e-11)
})

test_that("a hurdle interval (0, y2] is the base's mass exactly", {
  # it was (hu + tiny) - hu on the probability scale
  dp <- list(mu = 2, shape = exp(4), hu = 0.21)
  y2 <- c(0.5, 0.7, 1)
  got <- frmtmb:::row_lpdf(win_fam(hurdle_gamma()), c(0, 0, 0), c(0, 0, 0),
                           lapply(dp, rep, 3),
                           list(cens = rep(2, 3), cens_y2 = y2), NULL)
  ref <- log1p(-dp$hu) + stats::pgamma(y2, dp$shape, dp$shape / dp$mu,
                                       log.p = TRUE)
  expect_equal(got, ref, tolerance = 1e-12)
  # and on a count, where the base is truncated to above zero
  dq <- list(mu = 3, hu = 0.4)
  gq <- frmtmb:::row_lpdf(win_fam(hurdle_poisson()), c(0, 2), c(0, 2),
                          lapply(dq, rep, 2),
                          list(cens = c(2, 2), cens_y2 = c(1, 4)), NULL)
  p0 <- stats::dpois(0, 3)
  rq <- log1p(-dq$hu) - log1p(-p0) +
    log(stats::ppois(c(1, 4), 3) - stats::ppois(c(0, 2), 3))
  expect_equal(gq, rq, tolerance = 1e-12)
})

test_that("interval rows' AD gradient is the finite difference's", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  # at points away from the optimum, where a row's mass is about
  # exp(-35): there the probability-scale difference made the objective
  # rounding noise, and its gradient 0.2 to 13 off in relative terms
  cases <- list(
    skew_normal = list(fam = skew_normal(), y = c(1, 2.2, 2.9, 3.2),
                       w = 0.7, dp = list(mu = 0, sigma = 1, alpha = -10)),
    von_mises = list(fam = von_mises(), y = c(0.6, 1, 1.15), w = 0.4,
                     dp = list(mu = -2, kappa = 50)),
    hurdle_gamma = list(fam = hurdle_gamma(), y = c(0, 0, 1), w = 0.5,
                        dp = list(mu = 2, shape = exp(4), hu = 0.21)))
  for (nm in names(cases)) {
    k <- cases[[nm]]
    d <- data.frame(y = k$y, cc = 2, y2 = k$y + k$w)
    ob <- win_obj(k$fam, d, "y | cens(cc, y2)")
    p <- win_par(k$fam, k$dp)
    g <- as.numeric(ob$gr(p))
    fd <- numDeriv::grad(ob$fn, p, method = "Richardson")
    expect_equal(g, fd, tolerance = 1e-7, info = nm)
  }
})

test_that("trunc(lb = ) is divided by the log survivor at lb", {
  # with the mean far below the bound F(lb) rounds to one, and
  # log(1 - F(lb)) was -Inf: an infinite likelihood the optimizer could
  # step into
  set.seed(3)
  d <- data.frame(y = 10 + abs(stats::rnorm(50)))
  cases <- list(
    gaussian = list(fam = gaussian(), dp = list(mu = 0, sigma = 1),
                    ld = function(dp) stats::dnorm(d$y, dp$mu, dp$sigma,
                                                   log = TRUE),
                    ls = function(dp) stats::pnorm(10, dp$mu, dp$sigma,
                                                   lower.tail = FALSE,
                                                   log.p = TRUE)),
    student = list(fam = student(), dp = list(mu = 0, sigma = 1, nu = 4),
                   ld = function(dp) stats::dt(d$y - dp$mu, dp$nu,
                                               log = TRUE),
                   ls = function(dp) stats::pt(10 - dp$mu, dp$nu,
                                               lower.tail = FALSE,
                                               log.p = TRUE)),
    Gamma = list(fam = Gamma(link = "log"), dp = list(mu = 1, shape = 3),
                 ld = function(dp) stats::dgamma(d$y, dp$shape,
                                                 dp$shape / dp$mu,
                                                 log = TRUE),
                 ls = function(dp) stats::pgamma(10, dp$shape,
                                                 dp$shape / dp$mu,
                                                 lower.tail = FALSE,
                                                 log.p = TRUE)))
  for (nm in names(cases)) {
    k <- cases[[nm]]
    ob <- win_obj(k$fam, d, "y | trunc(lb = 10)")
    p <- win_par(k$fam, k$dp)
    expect_equal(ob$fn(p), -sum(k$ld(k$dp) - k$ls(k$dp)), tolerance = 1e-12,
                 info = nm)
    expect_true(all(is.finite(ob$gr(p))), info = nm)
  }
})

test_that("the weibull trunc() fuzz spec that stopped on a NaN gradient fits", {
  skip_on_cran()
  # fuzz spec 20446531: trunc(lb = 0.5) with a diag() term, one of the
  # three gated fuzz failures, was this normalizer
  plan <- fuzz_plan(seed = 20260901L, size = 300L)
  sp <- as.list(plan[plan$seed == 20446531L, , drop = FALSE][1, ])
  expect_identical(sp$family, "weibull")
  d <- fuzz_data(sp)
  fit <- eval(parse(text = fuzz_call_text(sp)))
  expect_equal(fit$opt$convergence, 0L)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("a custom family's cdf_log is read for every censored row", {
  # both tails on the log scale, and lcdf and lccdf derived from them
  my <- frmtmb_family(
    "my_normal", dpars = c("mu", "sigma"),
    links = list(mu = "identity", sigma = "log"),
    accepts_aterms = c("weights", "cens", "trunc"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dnorm(y, dpars[["mu"]], dpars[["sigma"]], log = TRUE)
    },
    cdf_log = function(q, dpars, aterms) {
      z <- (q - dpars[["mu"]]) / dpars[["sigma"]]
      list(lF = RTMB::pnorm(z, log.p = TRUE),
           lS = RTMB::pnorm(z, lower.tail = FALSE, log.p = TRUE))
    },
    init_dpars = list(mu = function(y, aterms) mean(y),
                      sigma = function(y, aterms) stats::sd(y)))
  expect_true(is.function(my$lcdf) && is.function(my$lccdf))
  set.seed(715)
  d <- data.frame(x = stats::rnorm(150))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(150)
  d$cc <- sample(c(-1, 0, 1), 150, TRUE)
  fit <- frm(bf(y | cens(cc) ~ x) + my, data = d)
  ref <- frm(y | cens(cc) ~ x, family = gaussian(), data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-10)
  expect_error(frmtmb_family("bad", dpars = "mu", links = list(mu = "log"),
                             lpdf = function(y, dpars, aterms) y,
                             cdf_log = 3),
               "cdf_log")
})

test_that("a summed CDF keeps its upper tail in a window", {
  # poisson at counts 1 to 30 and compois sum their masses, so log F near
  # 0 was accurate only in absolute terms: (30, 31] at mu 3 was -690.8
  dp <- list(mu = 3)
  w <- function(fam, dp, y, y2) {
    frmtmb:::row_lpdf(win_fam(fam), y, y, lapply(dp, rep, length(y)),
                      list(cens = rep(2, length(y)), cens_y2 = y2), NULL)
  }
  expect_equal(w(poisson(), dp, c(30, 20, 25), c(31, 25, 30)),
               c(stats::dpois(31, 3, log = TRUE),
                 log(sum(stats::dpois(21:25, 3))),
                 log(sum(stats::dpois(26:30, 3)))), tolerance = 1e-12)
  expect_equal(w(zero_inflated_poisson(), list(mu = 3, zi = 0.3), 30, 31),
               log(0.7) + stats::dpois(31, 3, log = TRUE), tolerance = 1e-12)
  expect_equal(w(hurdle_poisson(), list(mu = 3, hu = 0.3), 30, 31),
               log(0.7) + stats::dpois(31, 3, log = TRUE) -
                 log(-expm1(-3)), tolerance = 1e-12)
  lse <- function(v) max(v) + log(sum(exp(v - max(v))))
  lm <- RTMB::dcompois2(0:200, 3, 1.5, log = TRUE)
  expect_equal(w(brmsfamily("com_poisson"), list(mu = 3, nu = 1.5),
                 c(20, 31), c(25, 33)),
               c(lse(lm[22:26]), lse(lm[33:34])), tolerance = 1e-12)
})

test_that("trunc(lb = 31, ub = 60) on a poisson is exact away from its maximum", {
  # the window (30, 60] has its lower edge at the summed form's last
  # count; at an intercept of 0.5 the objective was +410624 for -35476
  set.seed(8)
  x <- stats::rnorm(4000)
  y <- stats::rpois(4000, exp(3.3 + 0.4 * x))
  keep <- y >= 31 & y <= 60
  d <- data.frame(x = x[keep], y = y[keep])
  ob <- frm(y | trunc(lb = 31, ub = 60) ~ x, family = poisson(), data = d,
            dry_run = "objective")$obj
  lse <- function(v) max(v) + log(sum(exp(v - max(v))))
  exact <- function(b) {
    mu <- exp(b[1] + b[2] * d$x)
    sum(stats::dpois(d$y, mu, log = TRUE) -
          vapply(mu, function(m) lse(stats::dpois(31:60, m, log = TRUE)), 0))
  }
  for (b0 in c(3.3, 2, 1, 0.5)) {
    expect_equal(-ob$fn(c(b0, 0.4)), exact(c(b0, 0.4)), tolerance = 1e-10,
                 info = format(b0))
  }
  fit <- frm(y | trunc(lb = 31, ub = 60) ~ x, family = poisson(), data = d)
  expect_equal(fit$opt$convergence, 0L)
  expect_equal(as.numeric(logLik(fit)), exact(fit$opt$par), tolerance = 1e-10)
})

test_that("compois interval rows have a smooth objective at large nu", {
  skip_if_not_installed("numDeriv")
  # with interval rows in the upper tail at nu = e^3 the central
  # differences were 3e5 and -3e7 where AD gave 1300
  set.seed(9)
  d <- data.frame(x = stats::rnorm(60))
  d$y <- stats::rpois(60, exp(1 + 0.3 * d$x))
  d$cc <- rep(c(0, 2), 30)
  d$y2 <- d$y + 2
  ob <- frm(bf(y | cens(cc, y2) ~ x, nu ~ 1),
            family = brmsfamily("com_poisson"), data = d,
            dry_run = "objective")$obj
  for (v in c(3, 4)) {
    p <- ob$par
    p[3] <- v
    expect_equal(as.numeric(ob$gr(p)),
                 numDeriv::grad(ob$fn, p, method = "Richardson"),
                 tolerance = 1e-6, info = format(v))
  }
})

test_that("simulate() fills every row of a valid truncated fit, from the window", {
  skip_on_cran()
  # the weibull fuzz spec 20446531 left 1 of 120 rows unfilled after 100
  # rejection passes; rows rejection cannot fill are drawn by inverting
  # the CDF
  plan <- fuzz_plan(seed = 20260901L, size = 300L)
  sp <- as.list(plan[plan$seed == 20446531L, , drop = FALSE][1, ])
  d <- fuzz_data(sp)
  fit <- eval(parse(text = fuzz_call_text(sp)))
  sims <- simulate(fit, nsim = 20, seed = 1)
  expect_false(anyNA(as.matrix(sims)))
  expect_true(all(as.matrix(sims) >= 0.5))
  # the row whose window keeps the least mass, drawn 4000 times, against
  # its truncated CDF: the KS distance under the 0.1 percent critical
  # value for that sample size
  fam <- frmtmb:::as_frmtmb_family(weibull())
  dp <- frmtmb:::eval_dpars(fit)[[1]]
  sc <- dp$mu / gamma(1 + 1 / dp$shape)
  keep <- stats::pweibull(0.5, dp$shape, sc, lower.tail = FALSE)
  i <- which.min(keep)
  m <- 4000L
  y <- frmtmb:::sim_response(fam, list(mu = rep(dp$mu[i], m),
                                       shape = rep(dp$shape[i], m)),
                             list(trunc_lb = rep(0.5, m)), m)
  Ft <- function(q) {
    1 - stats::pweibull(q, dp$shape[i], sc[i], lower.tail = FALSE) / keep[i]
  }
  ks <- max(abs(Ft(sort(y)) - seq_len(m) / m))
  expect_lt(ks, 1.95 / sqrt(m))
  expect_true(all(y >= 0.5))
})

test_that("inverse draws of a count stay in the window with its mean", {
  # a poisson at mu 2 truncated to [12, 15] keeps 1e-6 of its mass, so
  # rejection cannot fill it and every row is drawn by inversion
  fam <- frmtmb:::as_frmtmb_family(poisson())
  m <- 4000L
  set.seed(3)
  y <- frmtmb:::sim_response(fam, list(mu = rep(2, m)),
                             list(trunc_lb = rep(12, m), trunc_ub = rep(15, m)),
                             m)
  expect_true(all(y %in% 12:15))
  p <- stats::dpois(12:15, 2)
  p <- p / sum(p)
  expect_lt(abs(mean(y) - sum(12:15 * p)),
            4 * sqrt(sum((12:15)^2 * p) - sum(12:15 * p)^2) / sqrt(m))
})

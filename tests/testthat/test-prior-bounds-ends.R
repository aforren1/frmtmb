# Prior bounds at the ends of a density's support (backlog 5005, 5027)
# and on a centered intercept (5646). dev/priors2-bounds.R has the cases
# on 0.82.0 and on this build.

pbe_data <- function(seed = 1, n = 120) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n),
                  g = factor(rep(seq_len(n / 12), each = 12)))
  d$y <- 1 + d$z * 0.3 + stats::rnorm(n / 12)[d$g] * 0.3 +
    stats::rnorm(n, sd = 0.3)
  # the data put the slope of z near -0.5
  d$yneg <- d$y - 0.8 * d$z
  d
}

pbe_slope <- function(f) unname(fixef(f)["z", "Estimate"])

test_that("a bound where the density is zero is never the stop", {
  d <- pbe_data()
  for (dn in c("lognormal(0, 1)", "beta(2, 2)", "gamma(2, 2)")) {
    ub <- if (startsWith(dn, "beta")) 1 else NA
    # 0.82.0 stopped at z = 0 with the objective Inf and code 1. The box
    # run, which comes first, steps onto the end on its way, and its
    # "NA/NaN function evaluation" is about a point the refit left, so
    # it is not raised
    f0 <- frm(yneg ~ z, data = d,
              prior = set_prior(dn, class = "b", coef = "z", lb = 0,
                                ub = ub))
    expect_identical(f0$opt$convergence, 0L, label = dn)
    expect_true(is.finite(f0$opt$objective), label = dn)
    expect_gt(pbe_slope(f0), 0)
    # the mode the closed box finds once it starts inside the support,
    # a bound just above 0 that it does not reach
    f1 <- frm(yneg ~ z, data = d,
              prior = set_prior(dn, class = "b", coef = "z", lb = 1e-8,
                                ub = ub))
    expect_gt(pbe_slope(f1), 1e-6)
    expect_lt(abs(pbe_slope(f0) - pbe_slope(f1)) / pbe_slope(f1), 1e-5,
              label = dn)
    expect_lt(abs(f0$opt$objective - f1$opt$objective) /
                abs(f1$opt$objective), 1e-10, label = dn)
  }
})

test_that("an interior maximum is the one an unbounded fit finds", {
  d <- pbe_data()
  # a lognormal() slope with no bound is not boxed, and its maximum here
  # is inside; lb = 0 changes the chart, not the maximum
  free <- allow_warnings(
    frm(y ~ z, data = d,
        prior = set_prior("lognormal(0, 1)", class = "b", coef = "z")),
    "you have specified a lower bounded prior", require = "lower bounded")
  open <- frm(y ~ z, data = d,
              prior = set_prior("lognormal(0, 1)", class = "b", coef = "z",
                                lb = 0))
  expect_identical(open$opt$convergence, 0L)
  expect_lt(max(abs(open$opt$par - free$opt$par) /
                  pmax(1, abs(free$opt$par))), 1e-6)
})

test_that("an infinite-density end is refused only if the fit runs to it", {
  d <- pbe_data()
  # a gaussian slope the data put below 0: the likelihood is finite at
  # the bound, and the refit in the open chart runs to it (0.82.0
  # stopped there with the objective -Inf, code 1)
  expect_error(frm(yneg ~ z, data = d,
                   prior = set_prior("gamma(0.5, 1)", class = "b",
                                     coef = "z", lb = 0)),
               "z (gamma(0.5, 1), lb = 0). Toward that end",
               fixed = TRUE, class = "frmtmb_error")
  # beta(0.5, 2) is infinite at 0 and zero at 1: the refusal names the
  # end that is infinite, where it named ub = 1
  expect_error(frm(yneg ~ z, data = d,
                   prior = set_prior("beta(0.5, 2)", class = "b",
                                     coef = "z", lb = 0, ub = 1)),
               "z (beta(0.5, 2), lb = 0). Toward that end",
               fixed = TRUE, class = "frmtmb_error")
  # the guard absent, where the likelihood holds the fit off the end
  # (review of lane priors2, F1): a poisson("identity") slope, whose
  # likelihood is 0 at b = 0, has its maximum inside, where a bound
  # just off the end finds it; 0.82.0 failed with "NA/NaN gradient
  # evaluation"
  set.seed(11)
  dp <- data.frame(x = stats::runif(80, 0.5, 2))
  dp$y <- stats::rpois(80, 1.5 * dp$x)
  dp$y2 <- stats::rpois(80, 0.6 * dp$x)
  pois <- stats::poisson(link = "identity")
  for (cs in list(list(f = y ~ 0 + x, d = "gamma(0.5, 1)", ub = NA),
                  list(f = y2 ~ 0 + x, d = "beta(0.5, 2)", ub = 1))) {
    fe <- frm(cs$f, data = dp, family = pois,
              prior = set_prior(cs$d, class = "b", coef = "x", lb = 0,
                                ub = cs$ub))
    fi <- frm(cs$f, data = dp, family = pois,
              prior = set_prior(cs$d, class = "b", coef = "x", lb = 1e-10,
                                ub = cs$ub))
    expect_identical(fe$opt$convergence, 0L)
    expect_lt(abs(fe$opt$par - fi$opt$par) / abs(fi$opt$par), 1e-6,
              label = cs$d)
  }
  # a maximum inside the slack of the end is not a run to it (short
  # review of lane priors2, case H): an exposure of 1e7 per row puts it
  # at b* = (S - 1/2) / (sum(x) + 1) = 4.854e-10, and a first version
  # refused it on position alone
  set.seed(5)
  dh <- data.frame(x = stats::runif(80, 0.5, 2) * 1e7, y = 0L)
  dh$y[1L] <- 1L
  bstar <- 0.5 / (sum(dh$x) + 1)
  # its posterior vanishes at the end, so it is no local mode either
  fh <- allow_warnings(
    expect_no_warning(
      frm(y ~ 0 + x, data = dh, family = pois,
          prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                            lb = 0)),
      class = "frmtmb_prior_unbounded_local"),
    "Standard errors are not available")
  expect_identical(fh$opt$convergence, 0L)
  expect_lt(abs(unname(fh$opt$par) - bstar) / bstar, 1e-6)
  # a gaussian slope the data put inside: the likelihood is finite at
  # 0, so the posterior is unbounded there and the fit a local mode,
  # kept with a warning (0.82.0 ended at 0 with the objective -Inf)
  expect_warning(fy <-
    frm(y ~ z, data = d,
        prior = set_prior("gamma(0.5, 1)", class = "b", coef = "z",
                          lb = 0)),
    "z (gamma(0.5, 1), lb = 0)", fixed = TRUE,
    class = "frmtmb_prior_unbounded_local")
  expect_identical(fy$opt$convergence, 0L)
  expect_gt(pbe_slope(fy), 0.1)
  # the short review's case: y ~ x + z, the box run never near 0 (P3)
  set.seed(1)
  d3 <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150))
  d3$y <- 1 + 0.4 * d3$x - 0.3 * d3$z + stats::rnorm(150, sd = 0.7)
  expect_warning(f3 <-
    frm(y ~ x + z, data = d3,
        prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                          lb = 0)),
    "a local mode", fixed = TRUE, class = "frmtmb_prior_unbounded_local")
  expect_identical(f3$opt$convergence, 0L)
  # the guard absent: a density zero at the end, nothing said
  expect_no_warning(frm(y ~ x + z, data = d3,
                        prior = set_prior("gamma(2, 2)", class = "b",
                                          coef = "x", lb = 0)))
  # a standard deviation's bound at 0 is no box: its log is held, and
  # the prior's density on the log, times the Jacobian, vanishes at 0
  expect_s3_class(frm(y ~ z + (1 | g), data = d,
                      prior = set_prior("gamma(0.5, 1)", class = "sd",
                                        lb = 0)),
                  "frmtmb_fit")
  # the guard absent: a bound inside the support, and a density finite
  # at its end, keep their closed box and its maximum on the bound
  for (pr in list(set_prior("gamma(0.5, 1)", class = "b", coef = "z",
                            lb = 1e-8),
                  set_prior("exponential(1)", class = "b", coef = "z",
                            lb = 0),
                  set_prior("beta(1, 2)", class = "b", coef = "z", lb = 0,
                            ub = 1))) {
    # the fit's own rule says a bound holds the slope
    f <- allow_warnings(frm(yneg ~ z, data = d, prior = pr),
                        "Standard errors are not available",
                        require = "z: a bound holds it")
    expect_identical(f$opt$convergence, 0L)
  }
  # and the sampler is not refused: its posterior is proper
  expect_no_error(frmtmb:::resolve_prior_input(
    frm(y ~ z, data = d, dry_run = "objective"),
    set_prior("gamma(0.5, 1)", class = "b", coef = "z", lb = 0)))
})

test_that("a bound through a link lands inside the density's support", {
  bi <- frmtmb:::internal_bound
  # the link the fit reads a zi or a hu through
  lk <- frmtmb:::get_link("logit")
  lg <- list(scale = "natural", link = lk)
  ends <- seq(0.01, 0.99, by = 0.01)
  # the round trip through the link misses some of these ends
  rt <- lk$linkinv(lk$linkfun(ends))
  expect_gt(sum(rt != ends), 0)
  lo <- vapply(ends, bi, 0, pm = lg, which = "lb")
  hi <- vapply(ends, bi, 0, pm = lg, which = "ub")
  expect_true(all(lk$linkinv(lo) >= ends))
  expect_true(all(lk$linkinv(hi) <= ends))
  # by a few ulps, never more
  expect_lte(max(abs(lo - lk$linkfun(ends)) /
                   pmax(abs(lk$linkfun(ends)), 1e-300)),
             16 * .Machine$double.eps)
  sdp <- list(scale = "sd")
  ends2 <- seq(0.05, 5, by = 0.05)
  lo2 <- vapply(ends2, bi, 0, pm = sdp, which = "lb")
  hi2 <- vapply(ends2, bi, 0, pm = sdp, which = "ub")
  expect_true(all(exp(lo2) >= ends2) && all(exp(hi2) <= ends2))
  # near 1 the logit is flat in ulps of its argument, and a first
  # version that stepped one ulp at a time left 61 of 602 ends outside
  # (review of lane priors2, F3); each lands inside now
  ends3 <- c(1 - 10^-seq(2, 5, by = 0.01), 0.99981802991413904)
  lo3 <- vapply(ends3, bi, 0, pm = lg, which = "lb")
  hi3 <- vapply(ends3, bi, 0, pm = lg, which = "ub")
  expect_true(all(lk$linkinv(lo3) >= ends3))
  expect_true(all(lk$linkinv(hi3) <= ends3))
  # a bound no internal value maps inside is refused by name
  expect_error(frmtmb:::bound_inside(0, 2, stats::plogis, "lb"),
               "cannot be carried through the link", fixed = TRUE)
  # uniform(0, 0.3) on zi stopped at code 1, 0.22 short; now the
  # maximum is on the bound and the fit says so
  set.seed(1)
  dz <- data.frame(z = stats::rnorm(120))
  invisible(stats::rnorm(10))
  invisible(stats::rnorm(120))
  dz$cnt <- stats::rpois(120, 2) * stats::rbinom(120, 1, 0.5)
  f <- allow_warnings(frm(cnt ~ z, data = dz, family = zero_inflated_poisson(),
                          prior = set_prior("uniform(0, 0.3)", class = "zi")),
                      "Standard errors are not available",
                      require = "zi_(Intercept): a bound holds it")
  expect_identical(f$opt$convergence, 0L)
  zi <- stats::plogis(unname(f$estimates$betad[1L]))
  expect_lte(zi, 0.3)
  expect_lt((0.3 - zi) / 0.3, 1e-8)
  # a lower bound above the support maps to NaN, as one below it does;
  # 0.82.0 and a first version read both as "constrains nothing" and
  # fitted zi lb = 2 with no bound (short review of lane priors2)
  expect_error(frm(cnt ~ z, data = dz, family = zero_inflated_poisson(),
                   prior = set_prior("", class = "zi", lb = 2)),
               "lb = 2 is outside the support", fixed = TRUE,
               class = "frmtmb_error")
  # the guard absent: below the support the bound constrains nothing
  fz <- frm(cnt ~ z, data = dz, family = zero_inflated_poisson(),
            prior = set_prior("", class = "zi", lb = -1))
  f0 <- frm(cnt ~ z, data = dz, family = zero_inflated_poisson())
  expect_identical(fz$opt$objective, f0$opt$objective)
  # through the log: exp(log(0.35)) is below 0.35, and with the data's
  # sigma under it the fit stopped with "NA/NaN gradient evaluation";
  # 0.36 round-trips and was always fine
  set.seed(1)
  ds <- data.frame(z = stats::rnorm(120))
  ds$y <- 1 + 0.3 * ds$z + stats::rnorm(120, sd = 0.3)
  for (lb in c(0.35, 0.36)) {
    pr <- set_prior(sprintf("uniform(%.17g, 5)", lb), class = "sigma",
                    lb = lb, ub = 5)
    f <- allow_warnings(frm(y ~ z, data = ds, prior = pr),
                        "Standard errors are not available",
                        require = "sigma_(Intercept): a bound holds it")
    expect_identical(f$opt$convergence, 0L)
    s <- exp(unname(f$estimates$betad[1L]))
    expect_gte(s, lb)
    expect_lt((s - lb) / lb, 1e-8)
  }
})

test_that("a bound on a centered intercept is refused, as brms reads it", {
  d <- pbe_data()
  d$xs <- d$z + 10
  # brms 2.23.0 bounds its temporary intercept for centered predictors
  # (dev/priors2-brms.R, I1); 0.82.0 bounded the intercept at zero, and
  # the slope fell from 0.405 to 0.208
  for (pr in list(set_prior("", class = "Intercept", lb = 3),
                  set_prior("normal(0, 5)", class = "Intercept", ub = 9))) {
    expect_error(frm(y ~ xs, data = d, prior = pr),
                 "is a bound on the intercept at the means of the predictors",
                 fixed = TRUE, class = "frmtmb_error")
  }
  expect_error(frm(bf(y ~ 1, sigma ~ xs), data = d,
                   prior = set_prior("", class = "Intercept", dpar = "sigma",
                                     lb = -1)),
               "dpar = \"sigma\" is a bound", fixed = TRUE)
  # the guard absent: no predictor, a predictor scale() centered, and
  # frmtmb's spelling for the intercept at zero all take the bound
  fi <- allow_warnings(frm(y ~ 1, data = d,
                           prior = set_prior("", class = "Intercept",
                                             lb = 3)),
                       "Standard errors are not available")
  expect_equal(unname(fixef(fi)[1L, "Estimate"]), 3)
  d$zc <- as.numeric(scale(d$z))
  expect_s3_class(frm(y ~ zc, data = d,
                      prior = set_prior("", class = "Intercept", lb = -5)),
                  "frmtmb_fit")
  fb <- allow_warnings(frm(bf(y ~ xs, center = FALSE), data = d,
                           prior = set_prior("", class = "b",
                                             coef = "Intercept", lb = 3)),
                       "Standard errors are not available")
  expect_equal(unname(fixef(fb)["Intercept", "Estimate"]), 3)
})

test_that("a refusal at an infinite-density end needs the objective to fall", {
  probe <- frmtmb:::open_end_probe
  # a stand-in objective with no TMB state to restore
  fake <- function(f) list(env = new.env(), fn = f)
  bnd <- list(lower = 0, upper = Inf, inf_lower = TRUE, inf_upper = FALSE,
              end_label_lower = "b (gamma(0.5, 1), lb = 0)")
  at <- function(b, f) list(par = c(b = b), objective = f(b),
                            convergence = 0L, message = "")
  # the pole of the density outruns a likelihood finite at the end
  pole <- function(p) 0.5 * log(p[1L])
  r <- probe(fake(pole), pole, at(1e-9, pole), bnd, numeric(0))
  expect_identical(r$open_unbounded, "b (gamma(0.5, 1), lb = 0)")
  # a maximum inside the slack, case H's shape: the objective rises
  # toward the end, so nothing is said and nothing moves
  inner <- function(p) -0.5 * log(p[1L]) + 1e9 * p[1L]
  r <- probe(fake(inner), inner, at(5e-10, inner), bnd, numeric(0))
  expect_null(r$open_unbounded)
  expect_null(r$open_end_stop)
  expect_identical(unname(r$par), 5e-10)
  # an objective that falls toward the end but stays bounded: not a
  # proof, so the best point found is kept, at code 1, for a warning
  flat <- function(p) sqrt(p[1L])
  r <- probe(fake(flat), flat, at(1e-9, flat), bnd, numeric(0))
  expect_null(r$open_unbounded)
  expect_identical(r$open_end_stop, "b (gamma(0.5, 1), lb = 0)")
  expect_identical(r$convergence, 1L)
  expect_lt(unname(r$par), 1e-300)
  expect_warning(frmtmb:::check_map_unbounded(r, bnd),
                 "not shown to rise without bound", fixed = TRUE)
})

test_that("a stop on a bound at an infinite objective is not converged", {
  # the open chart keeps a zero-density end out of reach, so no fit of
  # this file stops there; the guard stays for a stop that does, a box
  # of a caller's own or a run that leaves the chart
  nf <- frmtmb:::nonfinite_not_converged
  r <- nf(list(par = c(a = 0, b = 1), objective = Inf, convergence = 0L,
               message = ""), lower = c(0, -Inf), upper = c(Inf, Inf))
  expect_identical(r$convergence, 1L)
  expect_match(r$message, "the objective is Inf")
  # the guard absent: a finite objective, or a stop inside the box
  r1 <- nf(list(par = c(a = 0, b = 1), objective = 3, convergence = 0L,
                message = ""), lower = c(0, -Inf), upper = c(Inf, Inf))
  expect_identical(r1$convergence, 0L)
  r2 <- nf(list(par = c(a = 0.5, b = 1), objective = Inf, convergence = 0L,
                message = ""), lower = c(0, -Inf), upper = c(Inf, Inf))
  expect_identical(r2$convergence, 0L)
})

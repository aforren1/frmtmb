# conditional_effects() on trunc() and se() terms whose variables the
# caller does not pin in `conditions`. brms's prepare_conditions() holds
# every variable of the model at its mean, those of the addition terms
# too, and the response with them: a trunc(lb = lo) bound is mean(lo),
# se(s) is mean(s), and trunc(lb = min(y) - 1) is mean(y) - 1, the
# expression on a grid whose y is mean(y). Up to 0.67.0 frmtmb refused
# each of these by name (ce_aterms()). brms 2.23.0 at frmtmb's
# estimates holds the same values and computes the same truncated mean
# (dev/fixes-ce-brms.R, dev/fixes-log/ce-brms-fix.txt); the block at
# the end of this file reruns that comparison.

cea_data <- local({
  set.seed(20261005)
  n <- 150
  d <- data.frame(x = stats::runif(n, -1, 1), lo = stats::runif(n, -1.5, 0),
                  s = stats::runif(n, 0.3, 1))
  d$y <- 1 + d$x + stats::rnorm(n)
  d <- d[d$y > d$lo & d$y < 10, ]
  d$y2 <- 1 + d$x + stats::rnorm(nrow(d), 0, d$s)
  d
})

cea_xs <- c(-1, 0, 1)

# sigma of a gaussian fit at the grid's x values
cea_sigma <- function(fit) {
  nd <- data.frame(x = cea_xs, lo = 0, s = 1)
  as.vector(frm_linpred(fit, newdata = nd, dpar = "sigma",
                        type = "response"))
}

cea_ce <- function(fit, method, ...) {
  suppressMessages(conditional_effects(
    fit, "x", method = method, int_conditions = list(x = cea_xs), ...))[[1]]
}

test_that("an unpinned bound or se is held at its mean, as brms holds it", {
  d <- cea_data
  ft <- frm(bf(y | trunc(lb = lo) ~ x), data = d)
  fe <- frm(bf(y | trunc(lb = min(y) - 1) ~ x), data = d)
  fs <- frm(bf(y2 | se(s) ~ x), data = d)
  for (m in c("posterior_epred", "posterior_predict")) {
    expect_identical(unique(cea_ce(ft, m, ndraws = 10)$lo), mean(d$lo))
    expect_identical(unique(cea_ce(fe, m, ndraws = 10)$y), mean(d$y))
    expect_identical(unique(cea_ce(fs, m, ndraws = 10)$s), mean(d$s))
  }
  # a pinned value is still the value used
  ce <- cea_ce(ft, "posterior_epred", conditions = list(lo = 0.25))
  expect_identical(unique(ce$lo), 0.25)
})

test_that("the expected value is the truncated mean at the held bound", {
  d <- cea_data
  for (lb_case in c("lo", "expr")) {
    fit <- if (lb_case == "lo") {
      frm(bf(y | trunc(lb = lo) ~ x), data = d)
    } else {
      frm(bf(y | trunc(lb = min(y) - 1) ~ x), data = d)
    }
    lb <- if (lb_case == "lo") mean(d$lo) else mean(d$y) - 1
    b <- fixef(fit)[, "Estimate"]
    mu <- b[["Intercept"]] + b[["x"]] * cea_xs
    sg <- cea_sigma(fit)
    a <- (lb - mu) / sg
    want <- mu + sg * stats::dnorm(a) / stats::pnorm(a, lower.tail = FALSE)
    expect_equal(cea_ce(fit, "posterior_epred")$estimate__, want,
                 label = lb_case)
  }
  fit <- frm(bf(y | trunc(ub = 10) ~ x), data = d)
  b <- fixef(fit)[, "Estimate"]
  mu <- b[["Intercept"]] + b[["x"]] * cea_xs
  sg <- cea_sigma(fit)
  bb <- (10 - mu) / sg
  expect_equal(cea_ce(fit, "posterior_epred")$estimate__,
               mu - sg * stats::dnorm(bb) / stats::pnorm(bb))
})

test_that("the predictive interval is the one of the held distribution", {
  d <- cea_data
  p <- c(0.025, 0.975)
  nsim <- 4000
  # the Monte Carlo error of a sample quantile, in units of the
  # distribution's own spread; the band must sit within 5 of it
  mc_ok <- function(got, q_exact, dens) {
    se <- sqrt(p * (1 - p) / nsim) / dens
    all(abs(got - q_exact) < 5 * se)
  }
  set.seed(7)
  ft <- frm(bf(y | trunc(lb = lo) ~ x), data = d)
  b <- fixef(ft)[, "Estimate"]
  mu <- b[["Intercept"]] + b[["x"]] * cea_xs
  sg <- cea_sigma(ft)
  ce <- cea_ce(ft, "posterior_predict", ndraws = nsim)
  pa <- stats::pnorm((mean(d$lo) - mu) / sg)
  for (k in seq_along(cea_xs)) {
    u <- pa[k] + p * (1 - pa[k])
    q <- mu[k] + sg[k] * stats::qnorm(u)
    dens <- stats::dnorm(q, mu[k], sg[k]) / (1 - pa[k])
    expect_true(mc_ok(c(ce$lower__[k], ce$upper__[k]), q, dens))
  }
  for (sig in c(FALSE, TRUE)) {
    fs <- if (sig) {
      frm(bf(y2 | se(s, sigma = TRUE) ~ x), data = d)
    } else {
      frm(bf(y2 | se(s) ~ x), data = d)
    }
    b <- fixef(fs)[, "Estimate"]
    mu <- b[["Intercept"]] + b[["x"]] * cea_xs
    sd_tot <- if (sig) sqrt(cea_sigma(fs)^2 + mean(d$s)^2) else
      rep(mean(d$s), 3)
    ce <- cea_ce(fs, "posterior_predict", ndraws = nsim)
    for (k in seq_along(cea_xs)) {
      q <- stats::qnorm(p, mu[k], sd_tot[k])
      expect_true(mc_ok(c(ce$lower__[k], ce$upper__[k]), q,
                        stats::dnorm(q, mu[k], sd_tot[k])), label = sig)
    }
  }
})

test_that("conditional_effects() on trunc() and se() equals brms's", {
  skip_unless_brms_fit()
  d <- cea_data
  cases <- list(
    list(bf(y | trunc(lb = lo) ~ x), brms::bf(y | trunc(lb = lo) ~ x)),
    list(bf(y | trunc(ub = 10) ~ x), brms::bf(y | trunc(ub = 10) ~ x)),
    list(bf(y | trunc(lb = min(y) - 1) ~ x),
         brms::bf(y | trunc(lb = min(y) - 1) ~ x)),
    list(bf(y2 | se(s) ~ x), brms::bf(y2 | se(s) ~ x)),
    list(bf(y2 | se(s, sigma = TRUE) ~ x),
         brms::bf(y2 | se(s, sigma = TRUE) ~ x)))
  for (cs in cases) {
    fit <- frm(cs[[1]], data = d)
    bb <- brms_fixed_fit(cs[[2]], gaussian(), d, fit, ndraws = 4)
    a <- brms::conditional_effects(bb, "x", method = "posterior_epred",
                                   int_conditions = list(x = cea_xs))[[1]]
    b <- cea_ce(fit, "posterior_epred")
    lab <- deparse1(cs[[1]]$formula)
    expect_exact_num(a$estimate__, b$estimate__, label = lab)
    for (v in intersect(c("lo", "s", "y", "y2"), names(b))) {
      expect_identical(unique(a[[v]]), unique(b[[v]]), label = v)
    }
  }
})

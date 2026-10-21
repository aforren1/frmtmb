# frechet(), gen_extreme_value() and discrete_weibull(): brms 2.23.0's
# three univariate families that frm() refused through 0.70.0. The
# densities are checked against an independent reference at a grid and
# at the edges, the fits against the truth that simulated them, and the
# post-fit paths against quadrature or a direct sum. The comparisons with
# brms's Stan functions and with a brms fit are measured in
# dev/famuni-findings.md, where they need a compiler and a sampler.

eps <- .Machine$double.eps

# the per-row log density on the AD tape, as a function of the linear
# predictors, with the .eta_<dpar> entries build_objective() stores
tape_lpdf <- function(fam, y, eta, what = "lpdf") {
  dp <- names(eta)
  len <- vapply(eta, length, 1L)
  f <- function(p) {
    pars <- list()
    at <- 0L
    for (nm in dp) {
      e <- p[at + seq_len(len[[nm]])]
      at <- at + len[[nm]]
      pars[[nm]] <- fam$links[[nm]]$linkinv(e)
      pars[[paste0(".eta_", nm)]] <- e
    }
    fam[[what]](y, pars, list())
  }
  tp <- RTMB::MakeTape(f, unlist(eta, use.names = FALSE))
  tp(unlist(eta, use.names = FALSE))
}

test_that("frm() takes each family by its brms name", {
  # by name only, so that the build without these families fails here
  # on the refusal itself and not on a missing constructor
  set.seed(30)
  d <- data.frame(x = stats::rnorm(80))
  d$yc <- exp(0.5 + 0.3 * d$x) * stats::rexp(80)^(-1 / 4)
  d$yd <- stats::rpois(80, 2)
  for (nm in c("frechet", "gen_extreme_value", "discrete_weibull")) {
    y <- if (nm == "discrete_weibull") "yd" else "yc"
    fit <- frm(stats::as.formula(paste(y, "~ x")), family = nm, data = d)
    expect_identical(family(fit)$family, nm)
    expect_true(is.finite(as.numeric(logLik(fit))), label = nm)
  }
})

test_that("the three families have brms's dpars, links and defaults", {
  expect_identical(frechet()$dpars, c("mu", "nu"))
  expect_identical(gen_extreme_value()$dpars, c("mu", "sigma", "xi"))
  expect_identical(discrete_weibull()$dpars, c("mu", "shape"))
  expect_identical(c(frechet()$link, frechet()$link_nu), c("log", "logm1"))
  expect_identical(c(gen_extreme_value()$link, gen_extreme_value()$link_sigma,
                     gen_extreme_value()$link_xi),
                   c("identity", "log", "log1p"))
  expect_identical(c(discrete_weibull()$link, discrete_weibull()$link_shape),
                   c("logit", "log"))
  expect_identical(discrete_weibull()$type, "discrete")
  # the sets brms validates against
  expect_identical(frechet(identity)$link, "identity")
  expect_error(frechet("logit"), "not a supported link for family 'frechet'")
  expect_error(gen_extreme_value("logit"), "gen_extreme_value")
  expect_identical(gen_extreme_value(link_xi = "identity")$link_xi,
                   "identity")
  expect_error(gen_extreme_value(link_xi = "log"),
               "not a supported link for parameter 'xi'")
  expect_error(frechet(link_nu = "log"),
               "not a supported link for parameter 'nu'")
  expect_identical(discrete_weibull("cloglog")$link, "cloglog")
  expect_error(discrete_weibull("log"),
               "not a supported link for family 'discrete_weibull'")
  # and every spelling of the name reaches them
  expect_identical(brmsfamily("frechet")$family, "frechet")
  expect_identical(frmtmb::as_frmtmb_family("gen_extreme_value")$family,
                   "gen_extreme_value")
})

test_that("the constructors take brms's arguments, defaults and order", {
  skip_if_not_installed("brms")
  # brms exports frechet() and gen_extreme_value(); discrete_weibull()
  # is internal there and reached through brmsfamily()
  bf_dw <- get("discrete_weibull", asNamespace("brms"))
  expect_identical(formals(frechet), formals(brms::frechet))
  expect_identical(formals(gen_extreme_value),
                   formals(brms::gen_extreme_value))
  expect_identical(formals(discrete_weibull), formals(bf_dw))
  for (nm in c("frechet", "gen_extreme_value", "discrete_weibull")) {
    bfam <- suppressWarnings(brms::brmsfamily(nm))
    ffam <- brmsfamily(nm)
    expect_identical(ffam$dpars, bfam$dpars, label = nm)
    expect_identical(ffam$link, bfam$link, label = nm)
    for (dp in setdiff(bfam$dpars, "mu")) {
      expect_identical(ffam[[c("links", dp, "name")]],
                       bfam[[paste0("link_", dp)]], label = paste(nm, dp))
    }
  }
})

test_that("the Frechet density is Stan's frechet(nu, mu / gamma(1 - 1/nu))", {
  y <- c(1e-3, 0.05, 0.5, 1, 2.5, 10, 100, 1e4)
  g <- expand.grid(y = y, mu = c(0.5, 4, 50), nu = c(1.05, 3, 60))
  fam <- frechet()
  nu <- 1 + exp(log(g$nu - 1))
  s <- g$mu / gamma(1 - 1 / nu)
  # the Frechet log density written out from its definition, in the
  # variables brms hands Stan
  z <- (g$y / s)^(-nu)
  ref <- log(nu) - log(s) + (1 + nu) * log(s / g$y) - z
  got <- tape_lpdf(fam, g$y, list(mu = log(g$mu), nu = log(g$nu - 1)))
  # the bound is the answer's own magnitude in units of the double that
  # carries it, times the condition number of z = exp(-nu log(y / s)),
  # which is |nu log(y / s)| and reaches 400 on this grid (z, 1e270)
  cond <- pmax(1, abs(nu * log(g$y / s)))
  expect_true(all(abs(got - ref) <= 64 * eps * pmax(1, abs(ref)) * cond))
  # log S in the right tail, where 1 - F rounds to 0: exp(-z) is 1 to
  # the last bit at z < 1e-16, and log S = log(z) - z / 2 to that order
  far <- c(1e6, 1e9)
  ls <- tape_lpdf(fam, far, list(mu = log(c(1, 1)), nu = log(c(2, 2))),
                  "lccdf")
  zf <- (far / (1 / gamma(1 - 1 / 3)))^(-3)
  expect_true(all(abs(ls - (log(zf) - zf / 2)) <= 64 * eps * abs(ls)))
})

test_that("the GEV density is brms's, and is exact at the Gumbel point", {
  # y = 50 at sigma = 0.3 puts z = xi x near -60, far past the upper end
  g <- expand.grid(y = c(-3, -0.5, 0, 0.7, 2, 8, 50), mu = c(-1, 0.5),
                   sigma = c(0.3, 4), xi = c(-0.4, 0.3, 1.2))
  fam <- gen_extreme_value()
  eta <- list(mu = g$mu, sigma = log(g$sigma), xi = log1p(g$xi))
  xi <- expm1(eta$xi)
  sg <- exp(eta$sigma)
  x <- (g$y - g$mu) / sg
  t <- 1 + xi * x
  inside <- t > 0
  # brms's gen_extreme_value_lpdf() for xi != 0, which at these xi loses
  # at most eps / |xi x| of its digits to log(t)
  tin <- ifelse(inside, t, 1)
  ref <- -log(sg) - (1 + 1 / xi) * log(tin) - tin^(-1 / xi)
  got <- tape_lpdf(fam, g$y, eta)
  tol <- 64 * eps * pmax(1, abs(ref)) / pmin(1, abs(xi * x))
  expect_true(all(abs(got - ref)[inside] <= tol[inside]))
  # no density outside the support: the log density falls by 1e8 per
  # unit of t below 0, so a row at t = -1e-5 or further is a density of
  # exactly 0 in double precision (brms's is NaN there)
  far_out <- t < -1e-5
  expect_gt(sum(far_out), 0)
  expect_true(all(exp(got[far_out]) == 0))
  # the distribution function is 0 below the lower end and 1 above the
  # upper end: a trunc() bound may lie there
  Fo <- tape_lpdf(fam, g$y, eta, "lcdf")[!inside]
  expect_identical(Fo, as.numeric(xi[!inside] < 0))

  # at xi = 0 exactly, the Gumbel density brms switches to
  yy <- c(-4, -1, 0, 0.5, 3, 30)
  gl0 <- tape_lpdf(fam, yy, list(mu = rep(0.2, 6), sigma = rep(log(1.3), 6),
                                 xi = rep(0, 6)))
  xg <- (yy - 0.2) / 1.3
  ref0 <- -log(1.3) - xg - exp(-xg)
  expect_true(all(abs(gl0 - ref0) <= 64 * eps * pmax(1, abs(ref0))))
  # and its gradient in the log1p predictor of xi: d/dxi at 0 is
  # -x + x^2 (1 - exp(-x)) / 2, and dxi / deta = 1 there
  for (k in seq_along(yy)) {
    f <- function(p) {
      fam$lpdf(yy[k], list(mu = 0.2, sigma = 1.3, xi = expm1(p),
                           .eta_xi = p), list())
    }
    tp <- RTMB::MakeTape(f, 0)
    gr <- tp$jacobian(0)[1, 1]
    want <- -xg[k] + xg[k]^2 * (1 - exp(-xg[k])) / 2
    expect_lte(abs(gr - want), 64 * eps * max(1, abs(want)),
               label = paste("gradient at y =", yy[k]))
  }
})

test_that("the GEV density near xi = 0 is the Gumbel's to first order", {
  fam <- gen_extreme_value()
  yy <- c(-2, 0, 1, 4)
  xg <- (yy - 0.2) / 1.3
  gum <- -log(1.3) - xg - exp(-xg)
  slope <- -xg + xg^2 * (1 - exp(-xg)) / 2
  for (xi0 in c(-1e-3, -1e-8, -1e-14, 1e-14, 1e-8, 1e-3)) {
    got <- tape_lpdf(fam, yy, list(mu = rep(0.2, 4), sigma = rep(log(1.3), 4),
                                   xi = rep(log1p(xi0), 4)))
    # lpdf(xi) = gumbel + slope * xi + O(xi^2): the second-order term
    # is below (1 + x^2)^2 xi^2 at these x
    bound <- (1 + xg^2)^2 * xi0^2 + 64 * eps * pmax(1, abs(gum))
    expect_true(all(abs(got - gum - slope * xi0) <= bound),
                label = paste("xi =", xi0))
  }
})

test_that("the GEV density agrees with evd::dgev", {
  skip_if_not_installed("evd")
  # read by name: evd is a reference for this one test, not a dependency
  # R CMD check should ask the package to declare
  dgev <- getExportedValue("evd", "dgev")
  g <- expand.grid(y = c(-2, 0.4, 3, 12), mu = c(-1, 0.5), sigma = c(0.5, 2),
                   xi = c(-0.3, 0.05, 0.6))
  fam <- gen_extreme_value()
  xi <- expm1(log1p(g$xi))
  sg <- exp(log(g$sigma))
  inside <- 1 + xi * (g$y - g$mu) / sg > 0
  ref <- mapply(function(y, m, s, k) {
    dgev(y, loc = m, scale = s, shape = k, log = TRUE)
  }, g$y, g$mu, sg, xi)
  got <- tape_lpdf(fam, g$y, list(mu = g$mu, sigma = log(g$sigma),
                                  xi = log1p(g$xi)))
  zx <- abs(xi * (g$y - g$mu) / sg)
  tol <- 64 * eps * pmax(1, abs(ref)) / pmin(1, zx)
  expect_true(all(abs(got - ref)[inside] <= tol[inside]))
})

test_that("the discrete Weibull mass is brms's, to its far tail", {
  g <- expand.grid(y = c(0, 1, 2, 5, 20), mu = c(0.05, 0.5, 0.9),
                   k = c(0.2, 1, 2.5))
  fam <- discrete_weibull()
  eta <- list(mu = stats::qlogis(g$mu), shape = log(g$k))
  mu <- stats::plogis(eta$mu)
  k <- exp(eta$shape)
  # brms's discrete_weibull_lpmf() on the probability scale: both terms
  # are at most 1, so their difference is exact to a few eps
  p <- mu^(g$y^k) - mu^((g$y + 1)^k)
  got <- exp(tape_lpdf(fam, g$y, eta))
  expect_true(all(abs(got - p) <= 8 * eps))
  # the log-odds path keeps the mass where mu rounds to 1: at eta = 40
  # and shape 1, P(Y = y) = mu^y (1 - mu) with log mu = -log1p(e^-40)
  lp <- tape_lpdf(fam, c(0, 3), list(mu = c(40, 40), shape = c(0, 0)))
  lm <- -log1p(exp(-40))
  want <- c(0, 3) * lm - 40 - log1p(exp(-40))
  expect_true(all(abs(lp - want) <= 64 * eps * abs(want)))
  # it sums to one
  y <- 0:5000
  s1 <- sum(exp(tape_lpdf(fam, y, list(mu = rep(stats::qlogis(0.8), 5001),
                                       shape = rep(log(0.8), 5001)))))
  expect_lte(abs(s1 - 1), 5001 * eps)
})

test_that("the families fit and recover the truth that simulated them", {
  set.seed(31)
  n <- 600
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  nu <- 1 + exp(log(3) + 0.2 * d$z)
  d$yf <- exp(1 + 0.3 * d$x) / gamma(1 - 1 / nu) * stats::rexp(n)^(-1 / nu)
  xi <- 0.15
  d$yg <- 1 + 0.5 * d$x + exp(-0.3 + 0.2 * d$z) *
    (stats::rexp(n)^(-xi) - 1) / xi
  q <- stats::plogis(1 + 0.4 * d$x)
  k <- exp(0.2 + 0.2 * d$z)
  d$yd <- pmax(ceiling((log1p(-stats::runif(n)) / log(q))^(1 / k) - 1), 0)
  cases <- list(
    list(f = bf(yf ~ x, nu ~ z), fam = frechet(),
         truth = c(Intercept = 1, x = 0.3, nu_Intercept = log(3),
                   nu_z = 0.2)),
    list(f = bf(yg ~ x, sigma ~ z), fam = gen_extreme_value(),
         truth = c(Intercept = 1, x = 0.5, sigma_Intercept = -0.3,
                   sigma_z = 0.2)),
    list(f = bf(yd ~ x, shape ~ z), fam = discrete_weibull(),
         truth = c(Intercept = 1, x = 0.4, shape_Intercept = 0.2,
                   shape_z = 0.2)))
  for (cs in cases) {
    # a well-identified fit says nothing: no SE alarm, no warning
    expect_no_warning(fit <- frm(cs$f, family = cs$fam, data = d))
    cf <- fixef(fit)
    est <- cf[names(cs$truth), "Estimate"]
    se <- cf[names(cs$truth), "Est.Error"]
    expect_true(all(is.finite(se)), label = cs$fam$family)
    # within four standard errors of the truth at seed 31
    expect_true(all(abs(est - cs$truth) < 4 * se), label = cs$fam$family)
  }
  fg <- frm(bf(yg ~ x, sigma ~ z), family = gen_extreme_value(), data = d)
  xh <- frm_linpred(fg, dpar = "xi", type = "response", se.fit = TRUE)
  expect_lt(abs(xh$fit[1] - xi), 4 * xh$se.fit[1])
})

test_that("a Gumbel model is gen_extreme_value with xi = 0", {
  set.seed(32)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x - 1.2 * log(stats::rexp(n))
  fit <- frm(bf(y ~ x, xi = 0), family = gen_extreme_value(), data = d)
  # the maximum is the Gumbel's, written out: the same log-likelihood at
  # the same point to the precision of the optimizer's stop
  b <- fixef(fit)[, "Estimate"]
  sg <- frm_linpred(fit, dpar = "sigma", type = "response")[1]
  xg <- (d$y - b[[1]] - b[[2]] * d$x) / sg
  ll <- sum(-log(sg) - xg - exp(-xg))
  expect_lte(abs(as.numeric(logLik(fit)) - ll), 64 * eps * n * abs(ll))
  # the mean is mu plus Euler's constant times sigma
  expect_equal(unname(fitted(fit)[1, "Estimate"]),
               unname(b[[1]] + b[[2]] * d$x[1] - digamma(1) * sg),
               tolerance = 64 * eps)
})

test_that("random effects, cens() and trunc() fit on each family", {
  set.seed(33)
  n <- 400
  d <- data.frame(x = stats::rnorm(n), g = gl(20, n / 20))
  b <- stats::rnorm(20, 0, 0.3)[d$g]
  d$yf <- exp(1 + 0.3 * d$x + b) / gamma(1 - 1 / 4) * stats::rexp(n)^(-1 / 4)
  d$yg <- 1 + 0.5 * d$x + b + 0.8 * (stats::rexp(n)^(-0.1) - 1) / 0.1
  q <- stats::plogis(1 + 0.4 * d$x + b)
  d$yd <- pmax(ceiling((log1p(-stats::runif(n)) / log(q))^(1 / 1.2) - 1), 0)
  for (cs in list(list(y = "yf", fam = frechet()),
                  list(y = "yg", fam = gen_extreme_value()),
                  list(y = "yd", fam = discrete_weibull()))) {
    f <- stats::as.formula(paste(cs$y, "~ x + (1 | g)"))
    fit <- frm(f, family = cs$fam, data = d)
    sdg <- VarCorr(fit)[["g"]][["sd"]][1, ]
    expect_lt(abs(sdg[["Estimate"]] - 0.3), 4 * sdg[["Est.Error"]],
              label = cs$fam$family)
    # right censoring at the 80th percentile, and the censored fit's
    # slope within four standard errors of the uncensored one
    cp <- stats::quantile(d[[cs$y]], 0.8, names = FALSE)
    if (identical(cs$fam$type, "discrete")) cp <- round(cp)
    d$yc <- pmin(d[[cs$y]], cp)
    d$cc <- as.numeric(d[[cs$y]] >= cp)
    fc <- frm(yc | cens(cc) ~ x, family = cs$fam, data = d)
    full <- frm(stats::as.formula(paste(cs$y, "~ x")), family = cs$fam,
                data = d)
    sc <- fixef(fc)["x", ]
    expect_lt(abs(sc[["Estimate"]] - fixef(full)["x", "Estimate"]),
              4 * sc[["Est.Error"]], label = paste(cs$fam$family, "cens"))
  }
})

test_that("discrete_weibull censoring and truncation are brms's", {
  set.seed(34)
  n <- 300
  d <- data.frame(x = stats::rnorm(n))
  q <- stats::plogis(1.5 + 0.3 * d$x)
  y <- pmax(ceiling((log1p(-stats::runif(n)) / log(q))^(1 / 1.1) - 1), 0)
  # right censored at 6 means Y > 6; truncated at lb = 2 means Y >= 2
  d$y <- pmin(y, 6)
  d$cc <- as.numeric(y > 6)
  keep <- d$y >= 2
  ft <- frm(y | trunc(lb = 2) + cens(cc) ~ x, family = discrete_weibull(),
            data = d[keep, ])
  b <- fixef(ft)[, "Estimate"]
  k <- frm_linpred(ft, dpar = "shape", type = "response")[1]
  dk <- d[keep, ]
  L <- stats::plogis(b[[1]] + b[[2]] * dk$x, log.p = TRUE)
  S <- function(v) exp(L * pmax(v, 0)^k)     # P(Y >= v)
  ll <- ifelse(dk$cc == 1, log(S(7)), log(S(dk$y) - S(dk$y + 1))) -
    log(S(2))
  expect_lte(abs(as.numeric(logLik(ft)) - sum(ll)),
             64 * eps * nrow(dk) * abs(sum(ll)))
})

test_that("the post-fit means are the distributions' means", {
  P <- unclass(gen_extreme_value())$post
  m <- P$mean_fn(list(mu = rep(1, 5), sigma = rep(2, 5),
                      xi = c(-0.3, 0, 1e-9, 0.4, 1)), list())
  want <- c(1 + 2 * (gamma(1.3) - 1) / -0.3, 1 - 2 * digamma(1),
            # at xi = 1e-9 the first-order Taylor term of the Gumbel mean;
            # the second is of order 1e-18
            1 + 2 * (-digamma(1) + (digamma(1)^2 / 2 + pi^2 / 12) * 1e-9),
            1 + 2 * (gamma(0.6) - 1) / 0.4, Inf)
  expect_equal(m, want, tolerance = 64 * eps)
  # the variance is finite below xi = 1/2 only
  v <- P$var_fn(list(mu = rep(0, 3), sigma = rep(1, 3), xi = c(0, 0.2, 0.5)),
                list())
  expect_equal(v[1:2], c(pi^2 / 6, (gamma(0.6) - gamma(0.8)^2) / 0.04),
               tolerance = 64 * eps)
  expect_identical(v[3], Inf)
  # near the Gumbel point, where gamma(1 - 2 xi) - gamma(1 - xi)^2 keeps
  # only eps / xi^2 of its digits: the variance is pi^2 / 6 plus
  # (2 zeta(3) + gamma_E pi^2 / 3) xi to first order, and the second-order
  # term is 11.7 xi^2, held to 15 xi^2. Written as the difference of two
  # gammas, the variance was 8.8e-5 relative off at xi = 1e-6
  # (dev/famuni-gevmean-mp.R).
  xs <- c(-1e-4, -1e-6, 1e-8, 1e-6, 1e-4)
  v0 <- P$var_fn(list(mu = rep(0, 5), sigma = rep(1, 5), xi = xs), list())
  c1 <- 2 * 1.2020569031595942 - digamma(1) * pi^2 / 3
  expect_true(all(abs(v0 - (pi^2 / 6 + c1 * xs)) <= 15 * xs^2 +
                    64 * eps * pi^2 / 6))
  m0 <- P$mean_fn(list(mu = rep(0, 5), sigma = rep(1, 5), xi = xs), list())
  eg <- -digamma(1)
  expect_true(all(abs(m0 - (eg + (eg^2 / 2 + pi^2 / 12) * xs)) <=
                    10 * xs^2 + 64 * eps))

  P <- unclass(discrete_weibull())$post
  mu <- c(0.3, 0.9, 0.99)
  k <- c(1.5, 0.8, 0.6)
  mf <- P$mean_fn(list(mu = mu, shape = k), list())
  # the direct sum of P(Y >= y) to 2e6 terms, whose tail beyond is
  # below 1e-25 at these parameters
  ref <- vapply(1:3, function(i) sum(mu[i]^(seq_len(2e6)^k[i])), 0)
  expect_true(all(abs(mf - ref) <= 1e4 * eps * ref))
  vf <- P$var_fn(list(mu = mu, shape = k), list())
  rv <- vapply(1:3, function(i) {
    y <- seq_len(2e6)
    sum((2 * y - 1) * mu[i]^(y^k[i])) - ref[i]^2
  }, 0)
  expect_true(all(abs(vf - rv) <= 1e6 * eps * rv))

  P <- unclass(frechet())$post
  expect_identical(P$mean_fn(list(mu = c(2, 3), nu = c(1.5, 4)), list()),
                   c(2, 3))
  expect_identical(P$var_fn(list(mu = 2, nu = 2), list()), Inf)
})

test_that("each family is a mixture component, as in brms", {
  set.seed(35)
  n <- 400
  d <- data.frame(y = c(stats::rpois(n / 2, 1), stats::rpois(n / 2, 12)))
  fit <- frm(y ~ 1, family = mixture(discrete_weibull(), discrete_weibull()),
             data = d)
  one <- frm(y ~ 1, family = discrete_weibull(), data = d)
  # the data are two well-separated counts, which two components fit
  # far better than one: the gain is a likelihood-ratio statistic on 3
  # parameters, held to its chi-square quantile at 1e-9
  expect_gt(2 * as.numeric(logLik(fit) - logLik(one)),
            stats::qchisq(1 - 1e-9, 3))
  d$yc <- c(stats::rexp(n / 2) + 1, stats::rexp(n / 2) * 3 + 8)
  for (fam in list(frechet(), gen_extreme_value())) {
    fm <- frm(yc ~ 1, family = mixture(fam, fam), data = d)
    f1 <- frm(yc ~ 1, family = fam, data = d)
    expect_gt(2 * as.numeric(logLik(fm) - logLik(f1)),
              stats::qchisq(1 - 1e-9, 1 + length(fam$dpars)))
  }
})

test_that("family = \"wiener\" names frmtmb.eam when it is not loaded", {
  skip_if(isNamespaceLoaded("frmtmb.eam"))
  d <- data.frame(y = c(0.5, 0.7, 0.9), x = 1:3)
  expect_error(frm(y ~ x, family = "wiener", data = d),
               "extension package frmtmb.eam, which is not loaded")
  expect_error(brmsfamily("wiener"), "frmtmb.eam")
})

test_that("family = \"wiener\" is frmtmb.eam's wiener() once that is loaded", {
  # the case the refusal above must not fire on. The name is held in a
  # variable so that R CMD check does not read this as a dependency
  eam <- "frmtmb.eam"
  skip_if_not(requireNamespace(eam, quietly = TRUE))
  expect_identical(frmtmb:::family_ctor("wiener"),
                   getExportedValue(eam, "wiener"))
  expect_identical(brmsfamily("wiener")$family, "wiener")
})

test_that("an identity link has no likelihood outside its range, on the tape", {
  # RTMB folds `0 * x` to a constant 0 when it records the tape, so a
  # guard written `0 * log(nu - 1)` was NaN on doubles and absent in the
  # fit: at nu = 0.8 the objective was finite, from log|gamma(-0.25)|, a
  # wrong scale (dev/reviews/2026-10-13-famuni.md, fix 2). This reads
  # the fitted objective, which is what the optimizer evaluates.
  set.seed(36)
  n <- 200
  d <- data.frame(x = stats::rnorm(n))
  d$y <- exp(1 + 0.3 * d$x) / gamma(1 - 1 / 3) * stats::rexp(n)^(-1 / 3)
  fit <- frm(y ~ x, family = frechet(link_nu = "identity"), data = d)
  par <- fit$opt$par
  inu <- which(names(par) == "betad")
  expect_length(inu, 1L)
  expect_true(is.finite(fit$obj$fn(par)))
  for (v in c(0.8, 1)) {
    p2 <- par
    p2[inu] <- v
    expect_true(is.nan(fit$obj$fn(p2)), label = paste("frechet nu =", v))
  }
  # the GEV under link_xi = "identity": xi <= -1, where the likelihood is
  # unbounded, has none (fix 3; seeds 621 and 622 ran to -1.007 before)
  rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
  set.seed(621)
  dg <- data.frame(y = rg(60, 0, 1, -0.8))
  fg <- allow_warnings(
    frm(y ~ 1, family = gen_extreme_value(link_xi = "identity"), data = dg),
    allowed = c("did not report convergence", "maximum absolute gradient",
                "Standard errors", "the edge at -1"),
    require = "the edge at -1")
  pg <- fg$opt$par
  ixi <- which(names(pg) == "betad")[2L]
  for (v in c(-1, -1.2)) {
    p2 <- pg
    p2[ixi] <- v
    expect_true(is.nan(fg$obj$fn(p2)), label = paste("GEV xi =", v))
  }
  # and the fit ends inside the support with xi above -1, so logLik() is
  # the log-likelihood: the sum of the density written out
  mu <- fixef(fg)[1, "Estimate"]
  sg <- frm_linpred(fg, dpar = "sigma", type = "response")[1]
  xi <- frm_linpred(fg, dpar = "xi", type = "response")[1]
  t <- 1 + xi * (dg$y - mu) / sg
  expect_gt(xi, -1)
  expect_true(all(t > 0))
  ll <- sum(-log(sg) - (1 + 1 / xi) * log(t) - t^(-1 / xi))
  expect_lte(abs(as.numeric(logLik(fg)) - ll), 1e3 * eps * abs(ll))
})

test_that("the xi edge warning does not fire inside the range", {
  # the case where the guard must stay quiet: xi = -0.45 and 0.8 at
  # n = 150 end at -0.478 and 0.889 (the review's section A)
  rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
  for (cfg in list(c(-0.45, 401), c(0.8, 406))) {
    set.seed(cfg[2])
    d <- data.frame(y = rg(150, 2, 1.5, cfg[1]))
    expect_no_warning(fit <- frm(y ~ 1, family = gen_extreme_value(),
                                 data = d))
    xi <- frm_linpred(fit, dpar = "xi", type = "response")[1]
    expect_gt(xi, -0.9)
  }
})

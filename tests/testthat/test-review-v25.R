# Heavy reference-validation file: full runs happen locally and in CI
# (NOT_CRAN=true). CRAN only needs to see the package work, not the
# validation program re-derived.
skip_on_cran()

# Confirmed defects from the v0.22-v0.25 code review.
#
# A1 is the one that changed answers: cens() and trunc() on the same
# response composed into a likelihood that was nobody's model. The rest
# are prediction- and sampling-surface bugs (non-conformable delta
# method under quadrature, inits outside Stan's bounds, silently
# dropped custom-family covariates, a discarded quadrature candidate, a
# raw LAPACK error, a decoder disagreeing with its own error message).

# --------------------------------------------------- A1 cens() x trunc()

# Draws from the population the fitted model claims: latent values
# truncated to [lb, ub] FIRST, then censored inside that window. This is
# also what simulate(censored = TRUE) generates.
sim_cens_trunc <- function(seed, n, lb, ub, cp, sigma = 0.9) {
  set.seed(seed)
  ys <- numeric(0)
  xs <- numeric(0)
  while (length(ys) < n) {
    xx <- stats::rnorm(2 * n)
    yy <- 0.4 + 0.5 * xx + stats::rnorm(2 * n, 0, sigma)
    ok <- yy > lb & yy < ub
    ys <- c(ys, yy[ok])
    xs <- c(xs, xx[ok])
  }
  data.frame(y = pmin(ys[seq_len(n)], cp), x = xs[seq_len(n)],
             cen = as.numeric(ys[seq_len(n)] > cp))
}

test_that("right censoring under trunc() is the windowed event probability", {
  dd <- sim_cens_trunc(2024, 400, -1, 3, 1.5)
  expect_gt(sum(dd$cen), 40)
  fit <- frm(bf(y | cens(cen) + trunc(lb = -1, ub = 3) ~ x) + gaussian(),
             data = dd)

  # P(y < Y <= ub | lb <= Y <= ub), not P(Y > y) / P(lb <= Y <= ub)
  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * dd$x
    s <- exp(p$ls)
    Z <- RTMB::pnorm((3 - mu) / s) - RTMB::pnorm((-1 - mu) / s)
    ll <- RTMB::dnorm(dd$y, mu, s, log = TRUE)
    ir <- which(dd$cen == 1)
    ll[ir] <- log(RTMB::pnorm((3 - mu[ir]) / s) -
                    RTMB::pnorm((dd$y[ir] - mu[ir]) / s))
    -sum(ll - log(Z))
  }
  est <- c(fixef(fit)$mu, exp = fixef(fit)$sigma[[1]])
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0), silent = TRUE)
  expect_lt(abs(obj$fn(unname(est)) + as.numeric(logLik(fit))), 1e-8)

  opt <- nlminb(unname(est), obj$fn, obj$gr)
  expect_lt(abs(as.numeric(logLik(fit)) - (-opt$objective)), 1e-6)
  expect_vector_equal(fixef(fit)$mu, opt$par[1:2], tol = 1e-4)

  # the pre-fix composition censored the UNtruncated variable and only
  # renormalized afterwards, which inflates the residual sd well past
  # the 0.9 the data were drawn with
  nll_old <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * dd$x
    s <- exp(p$ls)
    Z <- RTMB::pnorm((3 - mu) / s) - RTMB::pnorm((-1 - mu) / s)
    ll <- RTMB::dnorm(dd$y, mu, s, log = TRUE)
    ir <- which(dd$cen == 1)
    ll[ir] <- log(1 - RTMB::pnorm((dd$y[ir] - mu[ir]) / s))
    -sum(ll - log(Z))
  }
  o2 <- RTMB::MakeADFun(nll_old, list(b = c(0, 0), ls = 0), silent = TRUE)
  opt2 <- nlminb(unname(est), o2$fn, o2$gr)
  expect_gt(exp(opt2$par[3]) - exp(fixef(fit)$sigma[[1]]), 0.05)
  expect_lt(abs(exp(fixef(fit)$sigma[[1]]) - 0.9), 0.05)
})

test_that("left censoring under trunc() subtracts the truncation lower tail", {
  set.seed(2025)
  n <- 400
  ys <- numeric(0)
  xs <- numeric(0)
  while (length(ys) < n) {
    xx <- stats::rnorm(2 * n)
    yy <- 0.4 + 0.5 * xx + stats::rnorm(2 * n, 0, 0.9)
    ok <- yy > -1 & yy < 3
    ys <- c(ys, yy[ok])
    xs <- c(xs, xx[ok])
  }
  ys <- ys[seq_len(n)]
  xs <- xs[seq_len(n)]
  dd <- data.frame(y = pmax(ys, 0), x = xs, cen = -as.numeric(ys < 0))
  expect_gt(sum(dd$cen == -1), 40)
  fit <- frm(bf(y | cens(cen) + trunc(lb = -1, ub = 3) ~ x) + gaussian(),
             data = dd)

  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * dd$x
    s <- exp(p$ls)
    Flb <- RTMB::pnorm((-1 - mu) / s)
    Z <- RTMB::pnorm((3 - mu) / s) - Flb
    ll <- RTMB::dnorm(dd$y, mu, s, log = TRUE)
    il <- which(dd$cen == -1)
    ll[il] <- log(RTMB::pnorm((dd$y[il] - mu[il]) / s) - Flb[il])
    -sum(ll - log(Z))
  }
  est <- c(fixef(fit)$mu, fixef(fit)$sigma[[1]])
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0), silent = TRUE)
  expect_lt(abs(obj$fn(unname(est)) + as.numeric(logLik(fit))), 1e-8)
})

test_that("interval censoring under trunc() divides by the window mass", {
  set.seed(2026)
  n <- 300
  x <- stats::rnorm(n)
  ys <- 0.4 + 0.5 * x + stats::rnorm(n, 0, 0.9)
  ys <- pmin(pmax(ys, -0.9), 2.9)
  cen <- rep(0, n)
  cen[seq(1, n, by = 4)] <- 2
  y2 <- ys + 0.4
  y2[cen != 2] <- ys[cen != 2]
  dd <- data.frame(y = ys, x = x, cen = cen, y2 = pmin(y2, 2.95))
  fit <- frm(bf(y | cens(cen, y2) + trunc(lb = -1, ub = 3) ~ x) +
               gaussian(), data = dd)

  nll_ref <- function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    mu <- p$b[1] + p$b[2] * dd$x
    s <- exp(p$ls)
    Z <- RTMB::pnorm((3 - mu) / s) - RTMB::pnorm((-1 - mu) / s)
    ll <- RTMB::dnorm(dd$y, mu, s, log = TRUE)
    ii <- which(dd$cen == 2)
    ll[ii] <- log(RTMB::pnorm((dd$y2[ii] - mu[ii]) / s) -
                    RTMB::pnorm((dd$y[ii] - mu[ii]) / s))
    -sum(ll - log(Z))
  }
  est <- c(fixef(fit)$mu, fixef(fit)$sigma[[1]])
  obj <- RTMB::MakeADFun(nll_ref, list(b = c(0, 0), ls = 0), silent = TRUE)
  expect_lt(abs(obj$fn(unname(est)) + as.numeric(logLik(fit))), 1e-8)
})

test_that("the discrete censored-truncated form keeps the F(lb - 1) convention", {
  # cens() is refused for discrete families at the frame guard, so the
  # composed discrete branch is exercised on the objective directly
  set.seed(11)
  np <- 300
  xp <- stats::rnorm(np)
  yp <- stats::rpois(np, exp(0.3 + 0.4 * xp))
  keep <- yp >= 1 & yp <= 8
  dp <- data.frame(y = yp[keep], x = xp[keep])
  fp <- frm(bf(y | trunc(lb = 1, ub = 8) ~ x) + poisson(), data = dp)

  cen <- rep(0, nrow(dp))
  cen[seq(1, nrow(dp), by = 5)] <- 1
  cen[seq(2, nrow(dp), by = 7)] <- -1
  cen[dp$y >= 8 | dp$y <= 1] <- 0
  expect_gt(sum(cen == 1), 5)
  expect_gt(sum(cen == -1), 5)
  fr <- fp$frame
  fr$aterm_values$y$cens <- cen
  val <- frmtmb:::build_objective(fr)(fp$estimates)

  mu <- exp(fp$estimates$beta[1] + fp$estimates$beta[2] * dp$x)
  # inclusive lower bound: the window is P(1 <= Y <= 8) = F(8) - F(0),
  # and a left-censored count observes Y <= y, so F(y) needs no shift
  Z <- stats::ppois(8, mu) - stats::ppois(0, mu)
  ll <- stats::dpois(dp$y, mu, log = TRUE)
  ir <- cen == 1
  il <- cen == -1
  ll[ir] <- log(stats::ppois(8, mu[ir]) - stats::ppois(dp$y[ir], mu[ir]))
  ll[il] <- log(stats::ppois(dp$y[il], mu[il]) - stats::ppois(0, mu[il]))
  expect_lt(abs(val - (-sum(ll - log(Z)))), 1e-10)
})

test_that("OSA residuals stay calibrated on a cens() x trunc() fit", {
  dd <- sim_cens_trunc(77, 300, -1, 3, 1.5)
  fit <- frm(bf(y | cens(cen) + trunc(lb = -1, ub = 3) ~ x) + gaussian(),
             data = dd)
  r <- residuals(fit, type = "osa")
  obs <- dd$cen == 0
  expect_true(all(is.na(r[!obs])))
  expect_true(all(is.finite(r[obs])))

  # the OSA domain is the truncation window intersected with the
  # censoring window, [lb, cpoint]; an uncensored row's PIT is
  # (F(y) - F(lb)) / (F(cpoint) - F(lb)), which is what osa_cens_domain
  # and the corrected likelihood have to agree on
  dp <- frmtmb:::eval_dpars(fit)[["y"]]
  Flb <- stats::pnorm((-1 - dp$mu) / dp$sigma)
  Fcp <- stats::pnorm((1.5 - dp$mu) / dp$sigma)
  Fy <- stats::pnorm((dd$y - dp$mu) / dp$sigma)
  ref <- stats::qnorm((Fy - Flb) / (Fcp - Flb))
  expect_vector_equal(r[obs], ref[obs], tol = 1e-6)
  expect_gt(stats::ks.test(r[obs], "pnorm")$p.value, 0.01)
  expect_equal(frmtmb:::osa_cens_domain(fit$frame$aterm_values$y,
                                        fit$frame$y$y)$hi, 1.5)
})

test_that("the compat registry describes the composed cens x trunc likelihood", {
  rules <- frm_compat_rules()
  i <- which(rules$feature_a == "cens()" & rules$feature_b == "trunc()")
  expect_length(i, 1L)
  expect_identical(rules$status[i], "works")
  expect_match(rules$note[i], "truncated FIRST")
})

# ---------------------------------------------------- A2 deviance x se()

test_that("deviance residuals weight se() by the row's known variance", {
  set.seed(9)
  m <- 40
  sdv <- rep(c(0.1, 0.5), length.out = m)
  dd <- data.frame(x = stats::rnorm(m), sdy = sdv)
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(m, 0, sdv)

  fit <- frm(bf(y | se(sdy) ~ x) + gaussian(), data = dd)
  rd <- residuals(fit, type = "deviance")
  # se() alone maps sigma out at 1, so the prior weight is 1 / se^2 and
  # the deviance residual is the standardized residual
  expect_vector_equal(rd, (dd$y - fitted(fit)) / sdv, tol = 1e-8)
  # the defect: every row used to be scaled as if one dispersion
  # covered them all, so rd was exactly the raw residual
  raw <- dd$y - fitted(fit)
  scaling <- rd / raw
  expect_vector_equal(scaling, 1 / sdv, tol = 1e-8)
  expect_equal(max(scaling) / min(scaling), 5, tolerance = 1e-8)

  fitb <- frm(bf(y | se(sdy, sigma = TRUE) ~ x) + gaussian(), data = dd)
  dpb <- frmtmb:::eval_dpars(fitb)[["y"]]
  s_i <- sqrt(dpb$sigma^2 + sdv^2)
  expect_vector_equal(residuals(fitb, type = "deviance"),
                      (dd$y - fitted(fitb)) * dpb$sigma / s_i, tol = 1e-8)
})

test_that("the glm deviance agreement survives the se() weight", {
  set.seed(10)
  dd <- data.frame(x = stats::rnorm(60))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(60)
  fit <- frm(bf(y ~ x) + gaussian(), data = dd)
  ref <- stats::glm(y ~ x, data = dd)
  expect_vector_equal(residuals(fit, type = "deviance"),
                      unname(residuals(ref, type = "deviance")), tol = 1e-5)
})

# ------------------------------------------- B1 quadrature x se.fit

test_that("quadrature x predict(se.fit) reports modes-conditional SEs", {
  set.seed(4)
  ng <- 20
  nt <- 5
  n <- ng * nt
  dd <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                   x = stats::rnorm(n))
  dd$y <- stats::rpois(n, exp(0.3 + 0.4 * dd$x +
                                stats::rnorm(ng, 0, 0.5)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd,
             quadrature = TRUE)

  # the sdreport of a marginalized objective has no b rows, so the Z
  # block used to be cbind'd against an empty column-position vector
  expect_warning(p <- predict(fit, se.fit = TRUE),
                 "conditional on the conditional modes")
  expect_true(all(is.finite(p$se.fit)))
  expect_length(p$se.fit, n)
  expect_warning(pn <- predict(fit, newdata = dd, se.fit = TRUE),
                 "quadrature")
  expect_vector_equal(pn$se.fit, p$se.fit, tol = 1e-8)
  expect_warning(pr <- predict(fit, newdata = dd, se.fit = TRUE,
                               type = "response"), "quadrature")
  expect_true(all(is.finite(pr$se.fit)))

  # population-level prediction adds no b columns, so it is unaffected
  expect_silent(p0 <- predict(fit, se.fit = TRUE, re.form = NA))
  expect_true(all(is.finite(p0$se.fit)))
})

# ------------------------------------------------- B2 frm_sample bounds

test_that("mode_inits pulls every chain strictly inside the bounds", {
  lo <- c(0, -Inf, 2)
  hi <- c(Inf, 1, 2.5)
  set.seed(1)
  inits <- frmtmb:::mode_inits(c(-5, 5, 10), chains = 6, jitter = 2,
                               lower = lo, upper = hi)
  expect_length(inits, 6L)
  for (v in inits) {
    expect_true(all(v > lo))
    expect_true(all(v < hi))
  }
  # chain 1 is the mode anchor and is clamped like the rest
  expect_true(all(inits[[1]] > lo))
  # a zero jitter still clamps
  z <- frmtmb:::mode_inits(c(-5, 5, 10), chains = 2, jitter = 0,
                           lower = lo, upper = hi)
  expect_true(all(z[[1]] > lo))
  # unbounded sampling is untouched
  expect_equal(frmtmb:::mode_inits(c(1, 2), 1, 0)[[1]], c(1, 2))

  # a box narrower than the interior padding collapses to its midpoint
  expect_equal(frmtmb:::clamp_into_bounds(1e6, 1, 1 + 1e-9), 1 + 5e-10)
})

test_that("frm_sample warns when a bound excludes the ML mode", {
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")
  set.seed(55)
  n <- 120
  dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:12, 10)))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(12, 0, 0.4)[dd$g] +
    stats::rnorm(n)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  expect_lt(exp(fit$estimates$theta[[1]]), 1.5)

  # rstan's own message for an init at or outside a bound names neither
  # the parameter nor the bound, so frm_sample has to say it first
  seen <- character(0)
  ds <- withCallingHandlers(
    frm_sample(fit, chains = 2, iter = 300, refresh = 0, seed = 3,
               priors = set_prior("", class = "sd", lb = 1.5)),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("violates the requested bound", seen)))
  # the bound applies on the internal (log-sd) scale, and every draw
  # respects it: the chains started inside the box
  expect_true(all(as.matrix(ds)[, "theta_1"] >= log(1.5)))
})

# --------------------------------------------- B3 vint()/vreal() newdata

exposure_poisson <- function() {
  custom_family(
    "exposure_poisson",
    dpars = "mu",
    links = list(mu = "log"),
    lpdf = function(y, dpars, aterms) {
      lam <- dpars$mu * aterms$vreal1
      y * log(lam) - lam - lgamma(y + 1)
    },
    init_dpars = list(mu = function(y, aterms) mean(y) + 0.1),
    type = "discrete",
    post = list(
      mean_fn = function(dpars, aterms) dpars$mu * aterms$vreal1
    ),
    sim = function(dpars, aterms, n) {
      stats::rpois(n, dpars$mu * aterms$vreal1)
    }
  )
}

test_that("predict() requires the vreal() column instead of dropping it", {
  set.seed(31)
  n <- 200
  dd <- data.frame(x = stats::rnorm(n), expo = stats::runif(n, 0.5, 2))
  dd$y <- stats::rpois(n, exp(0.4 + 0.3 * dd$x) * dd$expo)
  fit <- frm(bf(y | vreal(expo) ~ x) + exposure_poisson(), data = dd)

  p <- predict(fit, newdata = dd, type = "response")
  expect_length(p, n)
  expect_vector_equal(p, unname(fitted(fit)), tol = 1e-8)

  # the payload used to be dropped, and the family's mean_fn then
  # returned a LENGTH-0 prediction with no message at all
  expect_error(predict(fit, newdata = data.frame(x = dd$x),
                       type = "response"),
               "vreal\\(expo\\).*no column expo")
  expect_error(frmtmb:::aterms_for_newdata(fit$spec$responses$y,
                                           data.frame(x = dd$x)),
               "no column expo")
})

# ---------------------------------------------------- B4 quad candidates

test_that("quad_fit keeps the lowest objective when none is stationary", {
  cand <- function(obj_value, stationary) {
    list(obj = TRUE, opt = list(objective = obj_value),
         stationary = stationary)
  }
  a1 <- cand(10, FALSE)
  a2 <- cand(3, FALSE)
  a3 <- cand(7, FALSE)
  best <- frmtmb:::quad_keep_best(NULL, a1)
  expect_identical(best, a1)
  best <- frmtmb:::quad_keep_best(best, a2)
  expect_identical(best$opt$objective, 3)
  best <- frmtmb:::quad_keep_best(best, a3)
  expect_identical(best$opt$objective, 3)

  # a stationary candidate wins whatever it costs, and a broken attempt
  # never displaces a good one
  expect_identical(frmtmb:::quad_keep_best(best, cand(50, TRUE))$stationary,
                   TRUE)
  expect_identical(frmtmb:::quad_keep_best(best, NULL), best)
  expect_identical(frmtmb:::quad_keep_best(best, list(retry = 1)), best)
})

# -------------------------------------------- B5 singular joint precision

test_that("predict(se.fit) degrades on a singular joint precision", {
  set.seed(21)
  n <- 120
  dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:12, 10)))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(12, 0, 0.4)[dd$g] +
    stats::rnorm(n)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = TRUE)
  p0 <- predict(fit, se.fit = TRUE)
  expect_true(all(is.finite(p0$se.fit)))

  # predict() inverted the joint precision by hand, so a singular one
  # threw a raw LAPACK message from inside the delta method
  Q <- fit$cache$sdr$jointPrecision
  sdr <- fit$cache$sdr
  sdr$jointPrecision <- Matrix::Matrix(1, nrow(Q), ncol(Q),
                                       dimnames = dimnames(Q))
  fit$cache$sdr <- sdr
  rm("Vjoint", envir = fit$cache)
  expect_warning(p <- predict(fit, se.fit = TRUE), "diagnose\\(\\)")
  expect_true(all(is.nan(p$se.fit)))
})

# ------------------------------------------------------ B6 cens() codes

test_that("cens() accepts the numeric codes its error message advertises", {
  expect_equal(frmtmb:::decode_cens(c("0", "1", "-1", "2")),
               c(0, 1, -1, 2))
  expect_equal(frmtmb:::decode_cens(factor(c("0", "1"))), c(0, 1))
  expect_equal(frmtmb:::decode_cens(c(" 1 ", "left")), c(1, -1))
  expect_error(frmtmb:::decode_cens(c("0", "zzz")), "cannot decode")

  set.seed(104)
  n <- 200
  x <- stats::rnorm(n)
  ystar <- 0.5 * x + stats::rnorm(n)
  y <- pmin(pmax(ystar, -1), 1.5)
  cen <- ifelse(ystar < -1, -1, ifelse(ystar > 1.5, 1, 0))
  dd <- data.frame(y = y, x = x, cnum = cen,
                   cstr = as.character(cen))
  f_num <- frm(bf(y | cens(cnum) ~ x) + gaussian(), data = dd)
  f_str <- frm(bf(y | cens(cstr) ~ x) + gaussian(), data = dd)
  expect_loglik_equal(f_num, f_str, tol = 1e-10)
})

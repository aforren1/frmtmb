## A centre-level frailty on a Royston-Parmar fit, measured rather than
## smoke tested.
##
## Item 2.5 of dev/extension-gaps-plan.md. The package used to assert
## that a `(1 | centre)` fit returned a finite log likelihood, which is
## not evidence of agreement with anything. What is asserted here is:
##
##   1. rstpm2 and frmtmb fit THE SAME model. The two write the spline
##      in different bases, so the check is a change of basis whose
##      residual is reported relative to the linear predictor's own
##      scale, and then an agreement in the estimates relative to the
##      run's own standard errors.
##   2. rstpm2's `theta` is the frailty VARIANCE, not its standard
##      deviation. The assertion fails if the two are swapped.
##   3. The gap between the two reported log likelihoods is the
##      INTEGRATION RULE and not the model: frmtmb's Laplace value at
##      its own optimum is compared with an exact adaptive-quadrature
##      integral written here in plain R, and the gap between the
##      packages is that difference.
##   4. A random effect on `gamma1` is a per-centre slope in log time,
##      it recovers, and it is anchored at t = 1, which is a unit and
##      not a fact.
##
## The sizes here are test sizes. The realistic design of the plan's
## "Realistic scale" table, 2000 subjects at 40 percent censoring, is
## measured in dev/frailty-findings.md over 200 replicates; every
## number below was read off a run whose script is named there.

# ------------------------------------------------------------ fixtures

## The Royston-Parmar basis, written here rather than reached for with
## `:::`, so that a disagreement with the package's own basis is
## visible instead of shared.
fr_basis <- function(knots, x) {
  nk <- length(knots)
  out <- cbind(rep(1, length(x)), x)
  if (nk > 2L) {
    kmin <- knots[1L]
    kmax <- knots[nk]
    for (j in seq_len(nk - 2L)) {
      kj <- knots[j + 1L]
      lam <- (kmax - kj) / (kmax - kmin)
      out <- cbind(out, pmax(x - kj, 0)^3 - lam * pmax(x - kmin, 0)^3 -
                     (1 - lam) * pmax(x - kmax, 0)^3)
    }
  }
  out
}

fr_dbasis <- function(knots, x) {
  nk <- length(knots)
  out <- cbind(rep(0, length(x)), rep(1, length(x)))
  if (nk > 2L) {
    kmin <- knots[1L]
    kmax <- knots[nk]
    for (j in seq_len(nk - 2L)) {
      kj <- knots[j + 1L]
      lam <- (kmax - kj) / (kmax - kmin)
      out <- cbind(out, 3 * pmax(x - kj, 0)^2 -
                     3 * lam * pmax(x - kmin, 0)^2 -
                     3 * (1 - lam) * pmax(x - kmax, 0)^2)
    }
  }
  out
}

## flexsurv's knot rule, which royston_parmar() also uses. Pinned here
## so that the two packages get the same knots by construction.
fr_knots <- function(logev, df) {
  bk <- range(logev)
  ik <- if (df > 1L) {
    unname(stats::quantile(logev, seq(0, 1, length.out = df + 1L)))[2:df]
  } else {
    numeric(0)
  }
  list(ik = ik, bk = bk, all = c(bk[1L], ik, bk[2L]))
}

## A Weibull proportional-hazards model with a shared log-normal
## frailty. That IS a Royston-Parmar model whose spline is linear in
## log time, so every truth is known: gamma0 = -shape log(scale),
## gamma1 = shape, every interior coefficient 0.
fr_sim <- function(seed, n = 800L, n_centre = 20L, sd_b = 0.5,
                   beta = 0.6, shape = 1.3, scale = 5, p_cens = 0.4) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  b <- stats::rnorm(n_centre, 0, sd_b)
  trt <- stats::rbinom(n, 1L, 0.5)
  eta <- beta * trt + b[centre]
  tt <- scale * (-log(stats::runif(n)) * exp(-eta))^(1 / shape)
  tau <- unname(stats::quantile(tt, 1 - p_cens))
  ev <- as.integer(tt <= tau)
  d <- data.frame(time = pmin(tt, tau), event = ev, censored = 1L - ev,
                  trt = trt, centre = factor(centre))
  attr(d, "truth") <- list(sd_b = sd_b, beta = beta, b = b,
                           gamma0 = -shape * log(scale), gamma1 = shape)
  d
}

## log H = gamma0 + beta trt + (gamma1 + u_c) log t: the centre carries
## its own Weibull SHAPE.
fr_sim_slope <- function(seed, n = 800L, n_centre = 20L, sd_u = 0.25,
                         gamma1 = 1.3, scale = 5, beta = 0.6,
                         p_cens = 0.4) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  u <- stats::rnorm(n_centre, 0, sd_u)
  trt <- stats::rbinom(n, 1L, 0.5)
  g0 <- -gamma1 * log(scale)
  tt <- exp((log(-log(stats::runif(n))) - g0 - beta * trt) /
              (gamma1 + u[centre]))
  tau <- unname(stats::quantile(tt, 1 - p_cens))
  ev <- as.integer(tt <= tau)
  d <- data.frame(time = pmin(tt, tau), event = ev, censored = 1L - ev,
                  trt = trt, centre = factor(centre))
  attr(d, "truth") <- list(sd_u = sd_u, gamma1 = gamma1, u = u,
                           shape = gamma1 + u, gamma0 = g0, beta = beta)
  d
}

## The marginal log likelihood of the shared log-normal frailty model,
## integrated per cluster by `stats::integrate()`. This is not any rule
## either package uses: frmtmb takes the Laplace approximation and
## rstpm2 takes adaptive Gauss-Hermite, and this is what both of them
## are approximations TO.
fr_exact_ll <- function(time, event, xb, cluster, knots, gam, sd_b) {
  x <- log(time)
  eta0 <- as.vector(fr_basis(knots, x) %*% gam) + xb
  lg <- log(as.vector(fr_dbasis(knots, x) %*% gam)) - x
  gfun <- function(b, idx) {
    e <- eta0[idx] + b
    sum(event[idx] * (e + lg[idx]) - exp(e)) +
      stats::dnorm(b, 0, sd_b, log = TRUE)
  }
  out <- 0
  for (idx in split(seq_along(time), cluster)) {
    op <- stats::optimize(function(b) -gfun(b, idx),
                          c(-12 * sd_b, 12 * sd_b), tol = 1e-12)
    bh <- op$minimum
    gh <- -op$objective
    h <- 1e-4
    curv <- (gfun(bh + h, idx) - 2 * gh + gfun(bh - h, idx)) / h^2
    sdl <- 1 / sqrt(max(-curv, 1e-8))
    iv <- stats::integrate(
      function(bs) vapply(bs, function(bb) exp(gfun(bb, idx) - gh),
                          numeric(1)),
      bh - 12 * sdl, bh + 12 * sdl, rel.tol = 1e-12,
      subdivisions = 2000L)
    out <- out + gh + log(iv$value)
  }
  out
}

fr_sd <- function(fit) sqrt(frmtmb::VarCorr(fit)[[1L]][1L, 1L])

fr_gam <- function(fit, df) {
  fx <- frmtmb::fixef(fit)
  c(unname(fx$mu[["(Intercept)"]]),
    vapply(paste0("gamma", seq_len(df)),
           function(p) unname(fx[[p]][["(Intercept)"]]), numeric(1)))
}

# ------------------------------------------------- the model, and whose

test_that("the frailty is the model rstpm2 fits, in another basis", {
  skip_if_not_installed("rstpm2")
  skip_if_not_installed("survival")
  skip_on_cran()
  d <- fr_sim(20260910L)
  kn <- fr_knots(log(d$time[d$event == 1L]), 2L)
  fit <- frmtmb::frm(
    frmtmb::bf(time | cens(censored) ~ trt + (1 | centre)),
    family = royston_parmar(knots = kn$ik, bknots = kn$bk), data = d,
    se = TRUE)
  # gsm() rather than its stpm2() wrapper: the wrapper rewrites its own
  # call to `gsm(...)` and evaluates it in the CALLER's frame, so it
  # needs rstpm2 attached, and a test file should not attach a
  # suggested package. gsm(penalised = FALSE) is what stpm2() calls.
  rst <- rstpm2::gsm(
    survival::Surv(time, event) ~ trt, data = d,
    smooth.formula = ~ rstpm2::nsx(log(time), knots = kn$ik,
                                   Boundary.knots = kn$bk),
    cluster = d$centre, RandDist = "LogN")

  ## THE CHANGE OF BASIS. rstpm2 writes the same natural cubic spline
  ## in the nsx basis; frmtmb writes it in Royston and Parmar's. Both
  ## span the natural cubic splines on these knots, so the map is
  ## linear and its residual is a check that the two spline SPACES are
  ## the same. If they were not, nothing below would mean anything.
  xg <- seq(kn$bk[1L] - 0.5, kn$bk[2L] + 0.5, length.out = 200L)
  cf <- stats::coef(rst)
  eta_r <- unname(cf[["(Intercept)"]]) +
    as.vector(rstpm2::nsx(xg, knots = kn$ik, Boundary.knots = kn$bk) %*%
                unname(cf[grep("nsx", names(cf))]))
  mp <- stats::lsfit(fr_basis(kn$all, xg), eta_r, intercept = FALSE)
  # relative to the linear predictor's own scale, never an absolute
  # number: measured 4.8e-15 at this seed
  expect_lt(max(abs(mp$residuals)), 1e-9 * max(abs(eta_r)))
  gam_r <- unname(mp$coefficients)

  ## THE ESTIMATES, each difference relative to the run's own standard
  ## error. Measured at this seed: 1.9e-04 of a standard error on the
  ## treatment coefficient.
  beta_f <- unname(frmtmb::fixef(fit)$mu[["trt"]])
  beta_r <- unname(cf[["trt"]])
  se_b <- summary(fit)[["coefficients"]][["mu"]]["trt", 2L]
  expect_lt(abs(beta_f - beta_r), 0.05 * se_b)
  se_br <- sqrt(rst@vcov["trt", "trt"])
  expect_lt(abs(se_b - se_br), 0.02 * se_b)

  ## THE FRAILTY SCALE. rstpm2 reports `logtheta` and theta is the
  ## VARIANCE of the log-normal frailty, so frmtmb's standard deviation
  ## matches its square root. Asserted as a RATIO so that reading theta
  ## as a standard deviation fails: measured 262 at this seed.
  th <- exp(unname(cf[["logtheta"]]))
  sd_f <- fr_sd(fit)
  expect_lt(abs(sd_f - sqrt(th)) * 50, abs(sd_f - th))
  # and on the log scale, against the interval the run itself reports
  cv <- frmtmb::confint_varcorr(fit)
  se_lsd <- (log(cv[["upr"]][1L]) - log(cv[["lwr"]][1L])) / (2 * 1.959964)
  expect_lt(abs(log(sd_f) - log(sqrt(th))), 0.1 * se_lsd)

  ## THE LOG LIKELIHOOD GAP IS THE INTEGRATION RULE. frmtmb integrates
  ## the frailty out by the Laplace approximation and rstpm2 by
  ## 9-node adaptive Gauss-Hermite, so the two report different numbers
  ## for the same model at almost the same point. The exact integral
  ## says which part of the gap is the rule: measured at this seed, the
  ## reported gap is -3.365e-02 and the Laplace offset alone is
  ## -3.361e-02, so 99.9 percent of it is the rule.
  ll_f <- as.numeric(stats::logLik(fit))
  # mle2 keeps the minimized objective, which is -logLik. logLik()
  # itself needs bbmle attached to dispatch, and the slot does not.
  ll_r <- -rst@min
  ex_f <- fr_exact_ll(d$time, d$event, beta_f * d$trt, d$centre, kn$all,
                      fr_gam(fit, 2L), sd_f)
  gap <- ll_f - ll_r
  expect_lt(abs(gap - (ll_f - ex_f)), 0.1 * abs(gap))
  ## and rstpm2's own rule is the closer one to the exact integral
  ex_r <- fr_exact_ll(d$time, d$event, beta_r * d$trt, d$centre, kn$all,
                      gam_r, sqrt(th))
  expect_lt(abs(ll_r - ex_r), 0.1 * abs(ll_f - ex_f))
})

test_that("the Laplace optimum is an optimum of the exact likelihood", {
  skip_on_cran()
  ## This needs no reference package at all. The Laplace approximation
  ## offsets the objective; the question a user cares about is whether
  ## it MOVES the maximizer, and that is answered by scoring the exact
  ## marginal likelihood at the fitted point and at points one quarter
  ## of a standard error away from it in each direction.
  d <- fr_sim(20260910L)
  kn <- fr_knots(log(d$time[d$event == 1L]), 2L)
  fit <- frmtmb::frm(
    frmtmb::bf(time | cens(censored) ~ trt + (1 | centre)),
    family = royston_parmar(knots = kn$ik, bknots = kn$bk), data = d,
    se = TRUE)
  beta <- unname(frmtmb::fixef(fit)$mu[["trt"]])
  sd_b <- fr_sd(fit)
  gam <- fr_gam(fit, 2L)
  se_b <- summary(fit)[["coefficients"]][["mu"]]["trt", 2L]
  cv <- frmtmb::confint_varcorr(fit)
  se_lsd <- (log(cv[["upr"]][1L]) - log(cv[["lwr"]][1L])) / (2 * 1.959964)
  at <- function(bb, ss) {
    fr_exact_ll(d$time, d$event, bb * d$trt, d$centre, kn$all, gam, ss)
  }
  best <- at(beta, sd_b)
  for (s in c(-0.25, 0.25)) {
    expect_lt(at(beta + s * se_b, sd_b), best)
    expect_lt(at(beta, sd_b * exp(s * se_lsd)), best)
  }
  ## and the Laplace offset is real rather than rounding: it is many
  ## multiples of the double precision resolution of a number this size
  expect_gt(abs(as.numeric(stats::logLik(fit)) - best),
            1e6 * .Machine$double.eps * abs(best))
})

# ------------------------------------------- the random effect on gamma1

test_that("a random effect on gamma1 is a per-centre shape, and it recovers", {
  skip_on_cran()
  ## gamma1 multiplies x = log t, so a centre deviation u makes that
  ## centre's cumulative hazard t^(gamma1 + u): a per-centre Weibull
  ## SHAPE, and a hazard ratio between two centres of t^(u - u'), which
  ## is time varying rather than proportional. It is NOT a frailty: a
  ## frailty multiplies the cumulative hazard by a constant.
  d <- fr_sim_slope(20260910L)
  tr <- attr(d, "truth")
  bk <- range(log(d$time[d$event == 1L]))
  fam <- royston_parmar(knots = numeric(0), bknots = bk)
  on <- suppressWarnings(frmtmb::frm(
    frmtmb::bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre)),
    family = fam, data = d, se = TRUE))
  off <- suppressWarnings(frmtmb::frm(
    frmtmb::bf(time | cens(censored) ~ trt), family = fam, data = d,
    se = TRUE))

  cv <- frmtmb::confint_varcorr(on)
  expect_true(cv[["lwr"]][1L] <= tr$sd_u && cv[["upr"]][1L] >= tr$sd_u)

  ## It has to beat the pooled fit at the thing it exists to estimate,
  ## which is the per-centre shape. Measured at this seed: 0.134
  ## against 0.236, a factor of 1.76.
  g1 <- unname(frmtmb::fixef(on)$gamma1[["(Intercept)"]])
  sh_hat <- g1 + as.numeric(frmtmb::ranef(on)[["centre"]])
  err_on <- mean(abs(sh_hat - tr$shape))
  err_off <- mean(abs(unname(frmtmb::fixef(off)$gamma1[["(Intercept)"]]) -
                        tr$shape))
  expect_lt(err_on, err_off)
  ## and the DEVIATIONS have to be what does it. Beating the pooled fit
  ## is not enough on its own: dropping this fit's own deviations and
  ## keeping its population gamma1 already beats the pooled fit at this
  ## seed, 0.234 against 0.236, so an inert block would pass the
  ## comparison above. Against its own population value the deviations
  ## are worth a factor of 1.74.
  expect_lt(err_on, 0.75 * mean(abs(g1 - tr$shape)))
  expect_gt(stats::cor(sh_hat, tr$shape), 0.5)

  ## Every centre's spline still has to increase in log time, or that
  ## centre has no hazard. The identity link on gamma1 does not hold it
  ## there, so the check is run rather than assumed.
  expect_equal(rp_floored(on, action = "report")$n_nonmonotone, 0L)
  expect_gt(min(sh_hat), 0)

  ## and that zero has to be a MEASURED zero. rp_floored() reads the
  ## fitted parameters, so a check that only ever saw the population
  ## gamma1 would report the same 0 on every fit. The floor was never
  ## reached by a fit in 75 tries (dev/frailty-findings.md), because at
  ## df = 1 the log density carries log(gamma1 + u) and that is a
  ## barrier for any centre with an event, so the case is constructed
  ## instead: push ONE centre's deviation past the population slope and
  ## the rows that come back must be that centre's own event rows.
  bad <- on
  j <- 3L
  # [["b"]] and not $b: names(fit$estimates) is beta, betad, b, theta,
  # and three of the four begin with "b". The exact match wins today,
  # and test-bracket-access.R exists because that is not a guarantee.
  bad$estimates[["b"]][j] <- -g1 - 0.2
  r <- rp_floored(bad, action = "report")
  expect_gt(r$n_nonmonotone, 0L)
  hit <- attr(r, "rows")[["nonmonotone"]]
  expect_setequal(as.character(d$centre[hit]), levels(d$centre)[j])
  expect_equal(r$n_nonmonotone,
               sum(d$centre == levels(d$centre)[j] & d$event == 1L))
  # the population value never left the monotone side, so nothing here
  # could have come from reading it
  expect_gt(g1, 0)
  expect_error(rp_floored(bad), "non-positive d\\(eta\\)/d\\(log t\\)")

  ## AND THE REACH OF THAT CHECK, asserted so the next reader is not
  ## told more than it does. rp_floored() tests `cens == 0` rows only
  ## (R/rp-check.R), so a group with NO events can carry a negative
  ## slope and be reported clean. Here that is shown by moving centre
  ## 3's own rows out of reach: with every one of them censored, the
  ## same broken deviation reports nothing at all.
  d2 <- d
  d2$censored[d2$centre == levels(d2$centre)[j]] <- 1L
  # exactly centre j's event rows moved and nothing else, so a 0 below
  # cannot come from having censored the whole dataset by accident
  expect_equal(sum(d2$censored) - sum(d$censored), r$n_nonmonotone)
  expect_gt(sum(d2$censored == 0L), 0L)
  bad2 <- bad
  bad2$frame[["aterm_values"]][["time"]][["cens"]] <-
    as.numeric(d2$censored)
  r2 <- rp_floored(bad2, action = "report")
  expect_equal(r2$n_nonmonotone, 0L)
  expect_silent(rp_floored(bad2))
  ## THIS ASSERTION PINS A GAP, NOT A GUARANTEE. It is here so that the
  ## reach of the check is written down rather than inferred, and it is
  ## meant to FAIL when the filed widening of `mono_rows` lands. Whoever
  ## lands it should flip these two lines to expect the count, not
  ## delete them.
  # dev/frailty-findings.md carries the fitted construction, where the
  # optimizer puts an all-censored centre below zero on its own.
})

test_that("a slope-only block on gamma1 is anchored at t = 1", {
  skip_on_cran()
  ## Rescale time by c and x = log t becomes x + log c, so
  ##   gamma0 + (gamma1 + u) x = (gamma0 - (gamma1 + u) log c) +
  ##                             (gamma1 + u) (x + log c)
  ## and the centre carries a random INTERCEPT of -u log c as well.
  ## A block on gamma1 alone therefore says every centre has the same
  ## cumulative hazard at t = 1 exactly, and t = 1 is a unit. A
  ## correlated block on both is the same model in any unit.
  ##
  ## The log likelihood of a density is not invariant to a rescale
  ## either: each EVENT row contributes log f and f scales by 1 / c, so
  ## the comparison is ll + n_event log c. Censored rows contribute
  ## log S, which is invariant.
  d <- fr_sim_slope(20260910L)
  nev <- sum(d$event)
  bk <- range(log(d$time[d$event == 1L]))
  mult <- 12
  at <- function(m, form) {
    dd <- d
    dd$time <- dd$time * m
    f <- suppressWarnings(frmtmb::frm(
      form, family = royston_parmar(knots = numeric(0),
                                    bknots = bk + log(m)),
      data = dd, se = TRUE))
    as.numeric(stats::logLik(f)) + nev * log(m)
  }
  slope <- frmtmb::bf(time | cens(censored) ~ trt,
                      gamma1 ~ (1 | centre))
  both <- frmtmb::bf(time | cens(censored) ~ trt + (1 | c | centre),
                     gamma1 ~ (1 | c | centre))
  d_slope <- abs(at(1, slope) - at(mult, slope))
  d_both <- abs(at(1, both) - at(mult, both))
  ## the paired block is invariant and the slope-only block is not, and
  ## the assertion is the RATIO between them rather than a size for
  ## either. Measured at the 2000-row design over a 365x rescale:
  ## 8.19 against 4.26e-08.
  expect_gt(d_slope, 100 * d_both)
})

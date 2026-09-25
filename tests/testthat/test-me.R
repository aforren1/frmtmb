# brms me(x, sdx): noise-free predictors with known measurement error
# (R/me.R). The reference for what the term MEANS is brms 2.23.0's Stan
# program (stan_Xme()); the reference for the fitted values is the
# closed-form marginal likelihood, which exists whenever the response is
# gaussian and every me() term enters linearly. There the latent values
# are jointly gaussian with the data, so the Laplace approximation is
# exact and logLik() must equal the multivariate-normal density.
#
# dev/me-findings.md reports the agreement with brms's own log density
# (dev/me-brms-lp.R).

me_data <- function(seed = 31, n = 90) {
  set.seed(seed)
  tx <- rnorm(n, 1, 0.8)
  tz <- 0.5 * tx + rnorm(n, 0, 0.7)
  sx <- runif(n, 0.2, 0.5)
  d <- data.frame(x = tx + rnorm(n, 0, sx), sx = sx,
                  z = tz + rnorm(n, 0, 0.3), sz = 0.3, w = rnorm(n))
  d$y <- 2 + 0.7 * tx - 0.4 * tz + 0.2 * d$w + rnorm(n, 0, 0.5)
  d$cnt <- rpois(n, exp(0.2 + 0.4 * tx))
  d
}

# The joint density of the data given the parameters, with the latent
# values integrated out analytically: each row is (y, x, z)
# multivariate normal. `b` holds the intercept, the slope on each latent
# value in `me_vars` order and the slope on w (0 when absent).
me_closed_ll <- function(d, b, bw, mu, sdv, R, sigma, me_vars) {
  M <- length(me_vars)
  S <- diag(sdv, M) %*% R %*% diag(sdv, M)
  out <- 0
  for (i in seq_len(nrow(d))) {
    noise <- vapply(me_vars, function(v) d[[paste0("s", v)]][i], 0)
    m <- c(b[1] + sum(b[-1] * mu) + bw * d$w[i], mu)
    Sb <- S %*% b[-1]
    V <- rbind(c(drop(t(b[-1]) %*% Sb) + sigma^2, Sb),
               cbind(Sb, S + diag(noise^2, M)))
    obs <- c(d$y[i], vapply(me_vars, function(v) d[[v]][i], 0))
    out <- out + mvtnorm::dmvnorm(obs, m, V, log = TRUE)
  }
  out
}

me_sigma <- function(fit) exp(fit$estimates[["betad"]][["sigma_(Intercept)"]])


test_that("one me() term: logLik is the closed-form marginal likelihood", {
  skip_if_not_installed("mvtnorm")
  d <- me_data()
  fit <- frm(bf(y ~ me(x, sx) + w) + gaussian(), data = d)
  fe <- fixef(fit)[, "Estimate"]
  ll <- me_closed_ll(d, b = fe[c("Intercept", "mexsx")], bw = fe[["w"]],
                     mu = fit$estimates$meanme[[1]],
                     sdv = exp(fit$estimates$logsdme[[1]]),
                     R = diag(1), sigma = me_sigma(fit), me_vars = "x")
  expect_lt(abs(as.numeric(logLik(fit)) - ll), 1e-10 * abs(ll))

  # and the maximum is the closed form's maximum: an independent
  # optimizer over the same six parameters lands on the same point
  nll <- function(p) {
    -me_closed_ll(d, b = p[1:2], bw = p[3], mu = p[4], sdv = exp(p[5]),
                  R = diag(1), sigma = exp(p[6]), me_vars = "x")
  }
  p0 <- c(fe[["Intercept"]], fe[["mexsx"]], fe[["w"]],
          fit$estimates$meanme[[1]], fit$estimates$logsdme[[1]],
          log(me_sigma(fit)))
  op <- stats::optim(p0 + 0.05, nll, method = "BFGS",
                     control = list(reltol = 1e-14, maxit = 2000))
  expect_lt(abs(op$value + as.numeric(logLik(fit))),
            1e-7 * abs(op$value))
  se <- fixef(fit)[c("Intercept", "mexsx", "w"), "Est.Error"]
  expect_lt(max(abs(op$par[1:3] - fe[c("Intercept", "mexsx", "w")]) / se),
            1e-3)

  # the correction the model exists for: the naive slope is attenuated
  naive <- stats::coef(stats::lm(y ~ x + w, data = d))[["x"]]
  expect_gt(fe[["mexsx"]], naive)
})

test_that("two me() terms are correlated by default, as in brms", {
  skip_if_not_installed("mvtnorm")
  d <- me_data(seed = 32)
  fit <- frm(bf(y ~ me(x, sx) + me(z, sz) + w) + gaussian(), data = d)
  fe <- fixef(fit)[, "Estimate"]
  hy <- summary(fit)$me
  expect_identical(rownames(hy), c("meanme_mex", "meanme_mez", "sdme_mex",
                                   "sdme_mez", "corme__mex__mez"))
  r <- hy["corme__mex__mez", "Estimate"]
  ll <- me_closed_ll(d, b = fe[c("Intercept", "mexsx", "mezsz")],
                     bw = fe[["w"]], mu = hy[c("meanme_mex", "meanme_mez"), 1],
                     sdv = hy[c("sdme_mex", "sdme_mez"), 1],
                     R = matrix(c(1, r, r, 1), 2), sigma = me_sigma(fit),
                     me_vars = c("x", "z"))
  expect_lt(abs(as.numeric(logLik(fit)) - ll), 1e-10 * abs(ll))

  # set_mecor(FALSE) fixes the correlation at zero: one parameter fewer,
  # and the closed form with R = I
  f0 <- frm(bf(y ~ me(x, sx) + me(z, sz) + w) + set_mecor(FALSE) +
              gaussian(), data = d)
  expect_null(f0$estimates[["thetame"]])
  expect_equal(attr(logLik(f0), "df") + 1, attr(logLik(fit), "df"))
  fe0 <- fixef(f0)[, "Estimate"]
  hy0 <- summary(f0)$me
  expect_false("corme__mex__mez" %in% rownames(hy0))
  ll0 <- me_closed_ll(d, b = fe0[c("Intercept", "mexsx", "mezsz")],
                      bw = fe0[["w"]],
                      mu = hy0[c("meanme_mex", "meanme_mez"), 1],
                      sdv = hy0[c("sdme_mex", "sdme_mez"), 1],
                      R = diag(2), sigma = me_sigma(f0),
                      me_vars = c("x", "z"))
  expect_lt(abs(as.numeric(logLik(f0)) - ll0), 1e-10 * abs(ll0))
  # the data were generated with correlated truths, which the
  # correlated model finds
  expect_gt(as.numeric(logLik(fit)), as.numeric(logLik(f0)))
})

test_that("me(gr = ) has one latent value per level", {
  skip_if_not_installed("mvtnorm")
  set.seed(33)
  ng <- 25
  g <- factor(sprintf("s%02d", rep(seq_len(ng), each = 4)))
  tx <- rnorm(ng, 0.5, 1)
  sxg <- runif(ng, 0.2, 0.6)
  xg <- tx + rnorm(ng, 0, sxg)
  d <- data.frame(g = g, xg = xg[g], sxg = sxg[g])
  d$y <- 1 + 0.8 * tx[g] + rnorm(nrow(d), 0, 0.4)
  fit <- frm(bf(y ~ me(xg, sxg, gr = g)) + gaussian(), data = d)
  expect_length(fit$estimates$miss, ng)
  fe <- fixef(fit)[, "Estimate"]
  b0 <- fe[["Intercept"]]
  b1 <- fe[["mexgsxggrEQg"]]
  mu <- fit$estimates$meanme[[1]]
  s <- exp(fit$estimates$logsdme[[1]])
  sig <- me_sigma(fit)
  ll <- 0
  for (j in seq_len(ng)) {
    rows <- which(as.integer(g) == j)
    k <- length(rows)
    V <- rbind(cbind(b1^2 * s^2 + diag(sig^2, k), b1 * s^2),
               c(rep(b1 * s^2, k), s^2 + sxg[j]^2))
    ll <- ll + mvtnorm::dmvnorm(c(d$y[rows], xg[j]),
                                c(rep(b0 + b1 * mu, k), mu), V, log = TRUE)
  }
  expect_lt(abs(as.numeric(logLik(fit)) - ll), 1e-10 * abs(ll))
  # brms's parameter names carry the level
  lab <- frmtmb:::brms_par_labels(fit)
  expect_true("Xme_mexg[s01]" %in% lab)

  # brms: x and sdx must be constant within a level
  d2 <- d
  d2$xg[2] <- d2$xg[2] + 1
  expect_error(frm(bf(y ~ me(xg, sxg, gr = g)) + gaussian(), data = d2),
               "should be unique for each group. Occured for level 's01'")
})

# The joint density brms's Stan program writes, at the fitted latent
# values, written out here from stan_Xme() and the likelihood rather
# than read off the tape. Where a me() term enters nonlinearly the
# Laplace approximation is no longer exact, so this is the check that
# the integrand itself is brms's.
me_joint_ll <- function(fit, d, ll_resp) {
  me <- fit$frame$me
  est <- fit$estimates
  miss <- est$miss
  out <- 0
  lat <- list()
  for (k in seq_along(me$terms)) {
    t <- me$terms[[k]]
    xl <- miss[t$idx]
    out <- out + sum(stats::dnorm(t$Xn, xl, t$noise, log = TRUE))
    lat[[k]] <- xl
  }
  for (grp in me$groups) {
    X <- do.call(cbind, lat[grp$K])
    mu <- est$meanme[grp$K]
    sdv <- exp(est$logsdme[grp$K])
    if (grp$cor) {
      C <- frmtmb:::us_chol_cor(est$thetame[grp$th_idx], length(grp$K))
      S <- diag(sdv) %*% C %*% diag(sdv)
      out <- out + sum(mvtnorm::dmvnorm(X, mu, S, log = TRUE))
    } else {
      for (j in seq_along(grp$K)) {
        out <- out + sum(stats::dnorm(X[, j], mu[j], sdv[j], log = TRUE))
      }
    }
  }
  out + ll_resp(lat)
}

test_that("interactions, dpars and non-gaussian families: the joint density", {
  skip_if_not_installed("mvtnorm")
  d <- me_data(seed = 34)
  joint <- function(fit) -fit$obj$env$f(fit$obj$env$last.par.best)

  f1 <- frm(bf(y ~ me(x, sx) * me(z, sz)) + gaussian(), data = d)
  fe <- fixef(f1)[, "Estimate"]
  expect_identical(names(fe),
                   c("Intercept", "mexsx", "mezsz", "mexsx:mezsz"))
  ref <- me_joint_ll(f1, d, function(l) {
    mu <- fe[["Intercept"]] + fe[["mexsx"]] * l[[1]] +
      fe[["mezsz"]] * l[[2]] + fe[["mexsx:mezsz"]] * l[[1]] * l[[2]]
    sum(stats::dnorm(d$y, mu, me_sigma(f1), log = TRUE))
  })
  expect_lt(abs(joint(f1) - ref), 1e-10 * abs(ref))

  f2 <- frm(bf(y ~ me(x, sx) * w) + gaussian(), data = d)
  fe <- fixef(f2)[, "Estimate"]
  expect_identical(names(fe), c("Intercept", "w", "mexsx", "mexsx:w"))
  ref <- me_joint_ll(f2, d, function(l) {
    mu <- fe[["Intercept"]] + fe[["w"]] * d$w + fe[["mexsx"]] * l[[1]] +
      fe[["mexsx:w"]] * l[[1]] * d$w
    sum(stats::dnorm(d$y, mu, me_sigma(f2), log = TRUE))
  })
  expect_lt(abs(joint(f2) - ref), 1e-10 * abs(ref))

  # two plain factors multiply; `:` between them is not R's sequence
  f5 <- frm(bf(y ~ me(x, sx):w:z) + gaussian(), data = d)
  fe <- fixef(f5)[, "Estimate"]
  expect_identical(names(fe), c("Intercept", "mexsx:w:z"))
  ref <- me_joint_ll(f5, d, function(l) {
    mu <- fe[["Intercept"]] + fe[["mexsx:w:z"]] * l[[1]] * d$w * d$z
    sum(stats::dnorm(d$y, mu, me_sigma(f5), log = TRUE))
  })
  expect_lt(abs(joint(f5) - ref), 1e-10 * abs(ref))

  f3 <- frm(bf(cnt ~ me(x, sx)) + poisson(), data = d)
  fe <- fixef(f3)[, "Estimate"]
  ref <- me_joint_ll(f3, d, function(l) {
    sum(stats::dpois(d$cnt, exp(fe[["Intercept"]] + fe[["mexsx"]] * l[[1]]),
                     log = TRUE))
  })
  expect_lt(abs(joint(f3) - ref), 1e-10 * abs(ref))

  # a dpar formula: the same latent values feed sigma
  f4 <- frm(bf(y ~ w, sigma ~ me(x, sx)) + gaussian(), data = d)
  b <- f4$estimates$beta
  bd <- f4$estimates$betad
  ref <- me_joint_ll(f4, d, function(l) {
    sum(stats::dnorm(d$y, b[["(Intercept)"]] + b[["w"]] * d$w,
                     exp(bd[["sigma_(Intercept)"]] + bd[["sigma_mexsx"]] *
                           l[[1]]), log = TRUE))
  })
  expect_lt(abs(joint(f4) - ref), 1e-10 * abs(ref))
  expect_true("bsp_sigma_mexsx" %in% variables(f4))
})

test_that("brms's names: bsp_, meanme_, sdme_, corme__ and Xme_", {
  d <- me_data(seed = 35)
  fit <- frm(bf(y ~ me(x, sx) * me(z, sz) + w) + gaussian(), data = d)
  v <- variables(fit)
  expect_true(all(c("bsp_mexsx", "bsp_mezsz", "bsp_mexsx:mezsz", "b_w",
                    "meanme_mex", "meanme_mez", "sdme_mex", "sdme_mez",
                    "corme__mex__mez") %in% v))
  # a hypothesis reaches the hyperparameters on brms's scale
  h <- hypothesis(fit, "sdme_mex > 0", class = NULL)
  expect_equal(h$hypothesis$Estimate,
               exp(fit$estimates$logsdme[["logsdme_mex"]]))
  lab <- frmtmb:::brms_par_labels(fit)
  expect_identical(lab[startsWith(lab, "Xme_mex[")],
                   paste0("Xme_mex[", seq_len(nrow(d)), "]"))
  # a transformed variable keeps brms's rename()
  f2 <- frm(bf(y ~ me(log(x + 5), sx)) + gaussian(), data = d)
  expect_true("bsp_melogxP5sx" %in% variables(f2))
  expect_true("sdme_melogxP5" %in% variables(f2))
})

test_that("fitted() uses the latent modes, newdata the observed values", {
  d <- me_data(seed = 36)
  fit <- frm(bf(y ~ me(x, sx) + w) + gaussian(), data = d)
  fe <- fixef(fit)[, "Estimate"]
  xl <- fit$estimates$miss[fit$frame$me$terms[[1]]$idx]
  expect_equal(unname(fitted(fit)[, "Estimate"]),
               unname(fe[["Intercept"]] + fe[["mexsx"]] * xl +
                        fe[["w"]] * d$w),
               tolerance = 1e-12)
  nd <- d[1:4, ]
  expect_equal(unname(fitted(fit, newdata = nd)[, "Estimate"]),
               unname(fe[["Intercept"]] + fe[["mexsx"]] * nd$x +
                        fe[["w"]] * nd$w),
               tolerance = 1e-12)
  # the latent modes shrink the noisy values toward meanme
  expect_lt(stats::var(xl), stats::var(d$x))
  expect_s3_class(simulate(fit, nsim = 2, seed = 1), "data.frame")
  nd$x[2] <- NA
  expect_error(fitted(fit, newdata = nd), "must supply complete numeric")
})

test_that("me() refusals name the term and the replacement", {
  d <- me_data(seed = 37)
  d$g <- factor(rep(1:9, 10))
  d$f <- factor(rep(c("a", "b", "c"), 30))
  d$o <- factor(sample(1:4, 90, TRUE), ordered = TRUE)
  # brms's own refusal (brms tests.brm.R:100)
  expect_error(frm(y ~ me(x, 2 * w) * me(x, w), data = d),
               "Variable 'x' is used in different calls to 'me'")
  expect_error(frm(y ~ me(x, sx) + me(x, sdx = sx), data = d),
               "Associated calls are: 'me(x, sx)', 'me(x, sdx = sx)'",
               fixed = TRUE)
  expect_error(frm(y ~ me(x), data = d), "Argument 'sdx' is missing")
  expect_error(frm(y ~ (me(x, sx) | g), data = d),
               "not supported inside a group-level term")
  expect_error(frm(y ~ I(me(x, sx)^2), data = d),
               "must be a term of its own or a factor")
  expect_error(frm(y ~ s(me(x, sx)), data = d),
               "must be a term of its own or a factor")
  expect_error(frm(y ~ mo(o):me(x, sx), data = d),
               "cannot share an interaction with mo[(][)] or mi[(][)]")
  expect_error(frm(y ~ me(x, sx):f, data = d),
               "me[(][)] interactions support numeric multipliers only")
  expect_error(frm(bf(y ~ a * me(x, sx), a ~ 1, nl = TRUE), data = d),
               "not supported in a nonlinear formula body")
  expect_error(frm(y ~ me(x, sx, gr = factor(g)), data = d),
               "takes the name of one grouping variable")
  d$s0 <- 0
  expect_error(frm(y ~ me(x, s0), data = d),
               "Measurement error should be positive")
  expect_error(frm(y ~ me(f, sx), data = d),
               "Noisy variables should be numeric")
  expect_error(frm(mvbf(bf(y ~ me(x, sx)) + set_mecor(FALSE),
                        bf(cnt ~ me(x, sx)) + set_mecor(TRUE)), data = d),
               "set_mecor[(][)] is TRUE on one response and FALSE")
  expect_error(set_mecor(NA), "must be TRUE or FALSE")

  # the fitting modes that cannot integrate the latent values
  expect_error(frm(y ~ me(x, sx) + (1 | g), data = d, quadrature = TRUE),
               "cannot be combined with me[(][)] or mi[(][)]")
  expect_error(frm(y ~ me(x, sx) + (1 | g), data = d, importance = 50),
               "cannot be combined with me[(][)] or mi[(][)]")
  fit <- frm(y ~ me(x, sx) + (1 | g), data = d)
  expect_error(vcov(fit, cluster = d$g),
               "does not support mi[(][)] / me[(][)] fits")
  # a likelihood with the measurements of x in it is not comparable
  # with one without them, even over the same rows
  f0 <- frm(y ~ x + (1 | g), data = d)
  expect_error(anova(fit, f0), "needs fits with the same me[(][)] terms")
  f1 <- frm(y ~ me(x, sx) + w + (1 | g), data = d)
  expect_s3_class(anova(fit, f1), "anova")
})

test_that("a multivariate model shares one latent vector per me() call", {
  d <- me_data(seed = 38)
  d$y2 <- 1 - 0.5 * d$x + rnorm(nrow(d), 0, 0.6)
  fit <- frm(bf(y ~ me(x, sx)) + bf(y2 ~ me(x, sx)) + gaussian(), data = d)
  expect_length(fit$estimates$miss, nrow(d))
  expect_true(all(c("bsp_y_mexsx", "bsp_y2_mexsx", "meanme_mex") %in%
                    variables(fit)))
  # set_mecor() on the multivariate formula reaches the spec
  f2 <- bf(y ~ me(x, sx) + me(z, sz)) + bf(y2 ~ me(x, sx))
  f2 <- f2 + set_mecor(FALSE)
  fit2 <- frm(f2 + gaussian(), data = d)
  expect_null(fit2$estimates[["thetame"]])
})

test_that("brms prior tables: the default rows drop, an edited one refuses", {
  skip_unless_brms()
  d <- me_data(seed = 39)
  bp <- with_brms_me(function() {
    brms::get_prior(y ~ me(x, sx) + me(z, sz), data = d)
  })
  expect_true(all(c("meanme", "sdme", "corme") %in% bp$class))
  # the table's Intercept and sigma rows are real priors and apply; its
  # me() rows are flat, lkj(1) included, and apply nothing
  f0 <- frm(bf(y ~ me(x, sx) + me(z, sz)) + gaussian(), data = d,
            prior = bp[!bp$class %in% c("meanme", "sdme", "corme"), ])
  f1 <- frm(bf(y ~ me(x, sx) + me(z, sz)) + gaussian(), data = d,
            prior = bp)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  bp$prior[bp$class == "sdme" & bp$coef == ""] <- "exponential(1)"
  expect_error(frm(bf(y ~ me(x, sx) + me(z, sz)) + gaussian(), data = d,
                   prior = bp), "brms's \"sdme\" is the SD of a me[(][)]")
  expect_error(set_prior("normal(0, 1)", class = "meanme"),
               "has no prior slot for it")
  # the me() coefficient is an ordinary class = "b" coefficient
  fp <- frm(bf(y ~ me(x, sx)) + gaussian(), data = d,
            prior = set_prior("normal(0, 0.01)", class = "b",
                              coef = "mexsx"))
  fm <- frm(bf(y ~ me(x, sx)) + gaussian(), data = d)
  expect_lt(abs(fixef(fp)["mexsx", 1]), abs(fixef(fm)["mexsx", 1]) / 5)
})

test_that("me() terms agree with brms's standata and prior names", {
  skip_unless_brms()
  set.seed(40)
  d <- data.frame(y = rnorm(40), x = rnorm(40), sx = runif(40, 0.1, 0.4),
                  z = rnorm(40), sz = 0.2, g = factor(rep(1:10, each = 4)))
  d$xg <- rep(rnorm(10), each = 4)
  d$sxg <- rep(runif(10, 0.1, 0.3), each = 4)
  fr <- frm(bf(y ~ me(x, sx) * me(z, sz) + me(xg, sxg, gr = g)) +
              gaussian(), data = d, dry_run = "frame")
  sd <- with_brms_me(function() {
    brms_standata(y ~ me(x, sx) * me(z, sz) + me(xg, sxg, gr = g), data = d)
  })
  ts <- fr$me$terms
  expect_identical(vapply(ts, `[[`, "", "coef"), c("mex", "mez", "mexg"))
  for (k in seq_along(ts)) {
    expect_equal(ts[[k]]$Xn, as.numeric(sd[[paste0("Xn_", k)]]))
    expect_equal(ts[[k]]$noise, as.numeric(sd[[paste0("noise_", k)]]))
  }
  expect_identical(fr$me$groups[[2]]$J, as.integer(sd$Jme_2))
  # the special-term coefficients, as brms's prior table spells them
  bp <- with_brms_me(function() {
    brms::get_prior(y ~ me(x, sx) * me(z, sz) + me(xg, sxg, gr = g),
                    data = d)
  })
  cn <- colnames(fr$linpreds[["y.mu"]]$X)
  expect_setequal(setdiff(cn, "(Intercept)"),
                  setdiff(bp$coef[bp$class == "b"], ""))
})

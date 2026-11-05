# A new row of a me() term has no latent value. brms 2.23.0 draws one per draw
# from normal(x, sdx) (prepare_predictions(), "sample new values of latent
# variables"). frm_linpred() here reads x and carries that draw's variance, and
# fitted() reads the mean over the draw. Seen to fail on 0.78.0 (rellib-r16): predict() on new data
# plugged x in, so its sd over draws was sigma alone, and fitted()'s Est.Error
# left the draw out (dev/predread-me-compare.R). brms's own
# numbers are in test-me-newdata-brms.R, and on draws in frmtmb.sample's
# test-readers-resp-me.R.

mn_data <- function() {
  set.seed(41)
  n <- 120
  g <- rep(1:12, each = 10)
  d <- data.frame(g = factor(g), z = stats::rnorm(n), x = stats::rnorm(n),
                  sx = 0.4, sw = 0.25, sg = 0.3)
  lat <- stats::rnorm(12)[g]
  d$xg <- lat + stats::rnorm(12, 0, 0.3)[g]
  d$w <- 0.6 * d$x + stats::rnorm(n)
  d$y <- 1 + 0.8 * d$x + 0.5 * d$w + 0.6 * lat + 0.3 * d$z +
    stats::rnorm(n, 0, 0.5)
  d
}
mn_grid <- function() {
  data.frame(x = c(-1, 0, 1), sx = c(0.2, 0.4, 0.8), w = 0.5, sw = 0.25,
             z = c(0, 1, -1), xg = c(0.5, 0.5, -1), sg = 0.3,
             g = factor(c("1", "1", "2"), levels = as.character(1:12)))
}
# one draw of the predictor, as a replicate reads it
mn_draw <- function(fit) {
  fit[["me_draw"]] <- new.env(parent = emptyenv())
  fit
}

test_that("a summary reads x and carries the draw's variance", {
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  nd <- mn_grid()
  b <- fixef(fit)
  bme <- b[grepl("^me", rownames(b)), "Estimate"]
  p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
  expect_equal(unname(p$fit), unname(b["Intercept", "Estimate"] +
                                       b["z", "Estimate"] * nd$z + bme * nd$x))
  lb <- frm_lp_basis(fit, newdata = nd, extra_cov = TRUE)
  expect_equal(lb$extra_var, (bme * nd$sx)^2)
  # rows without gr = are independent draws
  expect_equal(unname(as.matrix(lb$extra_cov)), diag((bme * nd$sx)^2))
  coef_var <- rowSums((lb$A %*% lb$V) * lb$A)
  expect_equal(unname(p$se.fit^2), unname(coef_var + (bme * nd$sx)^2))
  fi <- fitted(fit, newdata = nd)
  expect_equal(unname(fi[, "Est.Error"]), unname(p$se.fit))
  # on the fitted rows the latent values are the fit's own
  expect_equal(frm_lp_basis(fit)$extra_var, numeric(nrow(mn_data())))
})

test_that("rows at one level of me(gr = ) read one latent value", {
  fit <- frm(bf(y ~ me(xg, sg, gr = g) + z), data = mn_data(),
             family = gaussian())
  nd <- mn_grid()
  b <- fixef(fit)
  bme <- b[grepl("^me", rownames(b)), "Estimate"]
  lb <- frm_lp_basis(fit, newdata = nd, extra_cov = TRUE)
  Ec <- as.matrix(lb$extra_cov)
  vme <- (bme * 0.3)^2
  expect_equal(Ec[1, 2], vme)
  expect_equal(Ec[1, 3], 0)
  set.seed(1)
  e <- as.numeric(frm_linpred(mn_draw(fit), newdata = nd, dpar = "mu"))
  obs <- as.numeric(frm_linpred(fit, newdata = nd, dpar = "mu"))
  # one latent value: the rows moved together
  expect_equal(e[1] - obs[1], e[2] - obs[2])
  expect_false(isTRUE(all.equal(e[1] - obs[1], e[3] - obs[3])))
})

test_that("one draw reads one value, and draws differ", {
  fit <- frm(bf(y ~ me(x, sx) + me(x, sx):z), data = mn_data(),
             family = gaussian())
  nd <- transform(mn_grid(), z = 2)
  sh <- mn_draw(fit)
  set.seed(1)
  a <- frm_linpred(sh, newdata = nd)
  expect_identical(frm_linpred(sh, newdata = nd), a)
  expect_false(identical(frm_linpred(mn_draw(fit), newdata = nd), a))
  # both columns read one value of the call: |b1 + 2 b2| sdx, not their
  # independent sum. Over 4000 draws, 3 standard errors of a sd.
  b <- fixef(fit)
  rn <- rownames(b)[grepl("me", rownames(b))]
  expect_length(rn, 2L)
  e <- t(replicate(4000, frm_linpred(mn_draw(fit), newdata = nd)))
  want <- abs(b[rn[1L], "Estimate"] + 2 * b[rn[2L], "Estimate"]) * nd$sx
  ratio <- apply(e, 2, stats::sd) / want
  expect_lt(max(abs(ratio - 1)), 3 / sqrt(2 * 4000))
})

test_that("predict() and simulate() draw it per replicate", {
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  nd <- mn_grid()
  b <- fixef(fit)
  bme <- b[grepl("^me", rownames(b)), "Estimate"]
  sig <- sigma(fit)
  want <- sqrt(sig^2 + (bme * nd$sx)^2)
  n <- 6000L
  bound <- 4 / sqrt(2 * n)
  set.seed(2)
  pp <- predict(fit, newdata = nd, summary = FALSE, ndraws = n,
                propagate_error = FALSE)
  expect_lt(max(abs(apply(pp, 2, stats::sd) / want - 1)), bound)
  # the mean is the observed value's
  m <- frm_linpred(fit, newdata = nd)
  expect_lt(max(abs(colMeans(pp) - m) / (want / sqrt(n))), 4)
  s <- simulate(fit, nsim = n, seed = 3, newdata = nd)
  sa <- t(vapply(s, function(v) as.numeric(v), numeric(3L)))
  expect_lt(max(abs(apply(sa, 2, stats::sd) / want - 1)), bound)
  # the absent case: a model without me() is not touched, draw for draw
  f0 <- frm(bf(y ~ x + z), data = mn_data(), family = gaussian())
  a <- simulate(f0, nsim = 3, seed = 5, newdata = nd)
  sh <- f0
  sh[["me_draw"]] <- new.env(parent = emptyenv())
  expect_identical(a, simulate(sh, nsim = 3, seed = 5, newdata = nd))
  # on the fitted rows the replicate reads the fit's latent value
  set.seed(4)
  ins <- predict(fit, summary = FALSE, ndraws = 500L, propagate_error = FALSE)
  expect_lt(abs(stats::sd(ins[, 1L]) / sig - 1), 4 / sqrt(2 * 500))
})

test_that("conditional_effects(method = 'predict') draws it too", {
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  bme <- fixef(fit)[grepl("^me", rownames(fixef(fit))), "Estimate"]
  set.seed(5)
  ce <- conditional_effects(fit, "x", method = "predict", ndraws = 4000L,
                            resolution = 4L)[[1L]]
  want <- sqrt(sigma(fit)^2 + (bme * mean(mn_data()$sx))^2)
  expect_lt(max(abs(ce$se__ / want - 1)), 4 / sqrt(2 * 4000))
})

test_that("new data must carry sdx, positive, and one value per level", {
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  nd <- mn_grid()
  # brms 2.23.0 refuses each (dev/predread-me-brms-probe.R)
  expect_error(fitted(fit, newdata = nd[c("x", "z")]), "measurement SD")
  expect_error(predict(fit, newdata = nd[c("x", "z")]), "measurement SD")
  nd0 <- transform(nd, sx = 0)
  expect_error(fitted(fit, newdata = nd0), "should be positive")
  expect_error(predict(fit, newdata = nd0), "should be positive")
  fg <- frm(bf(y ~ me(xg, sg, gr = g) + z), data = mn_data(),
            family = gaussian())
  bad <- nd
  bad$xg[2L] <- 9
  expect_error(fitted(fg, newdata = bad), "unique for each group")
  expect_error(predict(fg, newdata = bad), "unique for each group")
  expect_error(fitted(fg, newdata = nd[c("xg", "sg", "z")]), "'g'")
  # a level the fit did not see, as a new level of a group-level term
  nw <- nd
  nw$g <- as.character(nw$g)
  nw$g[3L] <- "new"
  expect_error(fitted(fg, newdata = nw), "allow_new_levels",
               class = "frmtmb_new_levels")
  ok <- fitted(fg, newdata = nw, allow_new_levels = TRUE)
  expect_true(all(is.finite(ok[, "Est.Error"])))
  # a missing level is new too, each row its own
  na <- nd
  na$g[1:2] <- NA
  lb <- frm_lp_basis(fg, newdata = na, allow_new_levels = TRUE,
                     extra_cov = TRUE)
  expect_equal(as.matrix(lb$extra_cov)[1, 2], 0)
})

test_that("emmeans and the curve carry the variance, the derivative does not", {
  skip_if_not_installed("emmeans")
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  nd <- mn_grid()
  b <- fixef(fit)
  bme <- b[grepl("^me", rownames(b)), "Estimate"]
  em <- as.data.frame(emmeans::emmeans(fit, ~ z, at = list(z = c(0, 1)),
                                       cov.keep = character(0)))
  sx_bar <- mean(mn_data()$sx)
  x_bar <- mean(mn_data()$x)
  grid <- data.frame(x = x_bar, sx = sx_bar, z = c(0, 1))
  p <- frm_linpred(fit, newdata = grid, se.fit = TRUE)
  expect_equal(em$emmean, unname(p$fit))
  expect_equal(em$SE, unname(p$se.fit))
  cv <- frm_curve(fit, newdata = nd, simultaneous = FALSE)
  expect_equal(cv$.se,
               unname(frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit))
  d1 <- frm_curve_deriv(fit, var = "x", newdata = nd, simultaneous = FALSE)
  expect_equal(d1$.estimate, rep(unname(bme), 3L), tolerance = 1e-6)
  se_me <- unname(b[grepl("^me", rownames(b)), "Est.Error"])
  expect_equal(d1$.se, rep(se_me, 3L), tolerance = 1e-5)
})

# ---- the mean over the draw, and who reads it

mn_count <- function() {
  set.seed(7971)
  n <- 300
  d <- data.frame(z = stats::rnorm(n), lat = stats::rnorm(n), sx = 0.4)
  d$xo <- d$lat + stats::rnorm(n, 0, 0.4)
  eta <- 0.4 + 0.7 * d$lat
  d$cnt <- stats::rpois(n, exp(eta))
  d$pos <- exp(eta + stats::rnorm(n, 0, 0.5))
  d$hp <- ifelse(d$cnt == 0, 1L, d$cnt)
  d$hp[stats::runif(n) < 0.3] <- 0L
  d$bin <- stats::rbinom(n, 1, stats::plogis(eta))
  d$ord <- cut(eta + stats::rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
               labels = FALSE)
  d
}
mn_grid_c <- function() data.frame(xo = c(0, 1, 1), sx = c(0.8, 0.3, 1.2))

test_that("response_mean_gh() is the normal integral on the response scale", {
  fam <- as_frmtmb_family(poisson())
  eta <- c(-1, 0, 2)
  s2 <- c(0.2, 1, 4)
  got <- response_mean_gh(eta, s2, fam, list(mu = rep(1, 3)))
  expect_equal(got, exp(eta + s2 / 2), tolerance = 1e-10)
  # no variance is the value at the row
  expect_equal(response_mean_gh(eta, 0, fam, list(mu = rep(1, 3))), exp(eta))
  # a logit stays a probability
  fb <- as_frmtmb_family(bernoulli())
  gb <- response_mean_gh(c(1, 2), c(1, 9), fb, list(mu = c(0.5, 0.5)))
  expect_true(all(gb > 0.5 & gb < stats::plogis(c(1, 2))))
})

test_that("fitted() is the mean over the draw, and frm_linpred() stays at x", {
  d <- mn_count()
  nd <- mn_grid_c()
  fit <- frm(bf(cnt ~ me(xo, sx)), family = poisson(), data = d)
  b <- fixef(fit)
  b1 <- b[grepl("^me", rownames(b)), "Estimate"]
  eta <- b["Intercept", "Estimate"] + b1 * nd$xo
  fi <- fitted(fit, newdata = nd)
  # exp((b sdx)^2 / 2) times the value at x: the truth the mean must hit
  expect_equal(unname(fi[, "Estimate"]), exp(eta + (b1 * nd$sx)^2 / 2),
               tolerance = 1e-8)
  expect_equal(as.numeric(frm_linpred(fit, newdata = nd)), eta)
  expect_equal(as.numeric(frm_linpred(fit, newdata = nd, type = "response")),
               exp(eta))
  # the variance of a lognormal; the coefficients' part stays small beside it
  s2 <- (b1 * nd$sx)^2
  draw_sd <- exp(eta + s2 / 2) * sqrt(exp(s2) - 1)
  expect_true(all(fi[, "Est.Error"] >= draw_sd))
  expect_lt(max(fi[, "Est.Error"] / draw_sd - 1), 0.25)
  # the interval is the draw's quantiles, in the support, and holds the mean
  expect_true(all(fi[, "Q2.5"] > 0))
  expect_true(all(fi[, "Q2.5"] < fi[, "Estimate"] &
                    fi[, "Estimate"] < fi[, "Q97.5"]))
  cf <- frm_lp_basis(fit, newdata = nd)
  se2 <- rowSums((cf$A %*% cf$V) * cf$A)
  zq <- stats::qnorm(0.975)
  # the coefficients' share is read on the response scale, so on a log link it
  # sits near this form and not on it
  expect_equal(unname(fi[, "Q97.5"]), unname(exp(eta + zq * sqrt(se2 + s2))),
               tolerance = 0.005)
  expect_equal(unname(fi[, "Q2.5"]), unname(exp(eta - zq * sqrt(se2 + s2))),
               tolerance = 0.005)
  # against the draws themselves, 4000 of them
  set.seed(3)
  sims <- replicate(4000, {
    sh <- fit
    sh[["me_draw"]] <- new.env(parent = emptyenv())
    as.numeric(frm_linpred(sh, newdata = nd, type = "response"))
  })
  mcse <- apply(sims, 1, stats::sd) / sqrt(4000)
  expect_lt(max(abs(rowMeans(sims) - fi[, "Estimate"]) / mcse), 4.5)
})

test_that("a family whose mean is no one link inverse carries the draw", {
  d <- mn_count()
  nd <- mn_grid_c()
  specs <- list(list(f = bf(pos ~ me(xo, sx)), fam = lognormal()),
                list(f = bf(hp ~ me(xo, sx)), fam = hurdle_poisson()))
  for (m in specs) {
    fit <- frm(m$f, family = m$fam, data = d)
    fi <- fitted(fit, newdata = nd)
    # rows 2 and 3 differ in sdx alone
    expect_gt(fi[3L, "Est.Error"], 3 * fi[2L, "Est.Error"])
    expect_gt(fi[3L, "Estimate"], fi[2L, "Estimate"])
    expect_true(all(fi[, "Q2.5"] > 0))
    p <- frm_linpred(fit, newdata = nd, type = "response", se.fit = TRUE)
    expect_gt(p$se.fit[3L], 3 * p$se.fit[2L])
    band <- lapply(c(0.3, 1.2), function(s) {
      conditional_effects(fit, "xo", int_conditions = list(xo = 1),
                          conditions = list(sx = s))[[1L]]
    })
    expect_gt(band[[2L]]$upper__ - band[[2L]]$lower__,
              3 * (band[[1L]]$upper__ - band[[1L]]$lower__))
    expect_gt(band[[2L]]$lower__, 0)
    # the draws say the same
    set.seed(5)
    sims <- replicate(2000, {
      sh <- fit
      sh[["me_draw"]] <- new.env(parent = emptyenv())
      as.numeric(frm_linpred(sh, newdata = nd, type = "response"))
    })
    expect_lt(max(abs(rowMeans(sims) - fi[, "Estimate"]) /
                    (apply(sims, 1, stats::sd) / sqrt(2000))), 4.5)
    r <- fi[, "Est.Error"] / apply(sims, 1, stats::sd)
    expect_true(all(r > 0.8 & r < 1.3))
  }
})

test_that("a category distribution integrates, and its rows still sum to one", {
  d <- mn_count()
  nd <- mn_grid_c()
  fit <- frm(bf(ord ~ me(xo, sx)), family = cumulative(), data = d)
  fi <- fitted(fit, newdata = nd)
  expect_equal(unname(rowSums(fi[, "Estimate", ])), rep(1, 3))
  set.seed(5)
  sims <- replicate(2000, {
    sh <- fit
    sh[["me_draw"]] <- new.env(parent = emptyenv())
    frm_linpred(sh, newdata = nd, type = "response")[3L, ]
  })
  mcse <- apply(sims, 1, stats::sd) / sqrt(2000)
  expect_lt(max(abs(rowMeans(sims) - fi[3L, "Estimate", ]) / mcse), 4.5)
  expect_true(all(abs(fi[3L, "Est.Error", ] /
                        apply(sims, 1, stats::sd) - 1) < 0.25))
})

test_that("a product of two calls keeps its whole variance", {
  set.seed(4)
  d <- data.frame(x = stats::rnorm(150), w = stats::rnorm(150), sx = 0.3,
                  sw = 0.3)
  d$y <- 0.9 * d$x * d$w + stats::rnorm(150, 0, 0.5)
  fit <- frm(bf(y ~ me(x, sx):me(w, sw)), data = d, family = gaussian())
  nd <- data.frame(x = c(1, 0), w = c(0.5, 0), sx = 1.2, sw = 0.4)
  b1 <- fixef(fit)[2L, "Estimate"]
  want <- b1^2 * ((nd$x^2 + 1.2^2) * (nd$w^2 + 0.4^2) - nd$x^2 * nd$w^2)
  expect_equal(frm_lp_basis(fit, newdata = nd)$extra_var, want)
  set.seed(2)
  e <- t(replicate(4000, {
    sh <- fit
    sh[["me_draw"]] <- new.env(parent = emptyenv())
    as.numeric(frm_linpred(sh, newdata = nd))
  }))
  expect_lt(max(abs(apply(e, 2, stats::sd) / sqrt(want) - 1)), 0.06)
})

test_that("a display of me(gr = ) holds the group at NA, and answers", {
  fit <- frm(bf(y ~ me(xg, sg, gr = g) + z), data = mn_data(),
             family = gaussian())
  ce <- conditional_effects(fit, "xg", resolution = 4L)[[1L]]
  expect_true(all(is.na(ce$g)))
  expect_true(all(is.finite(ce$estimate__)))
  expect_setequal(names(conditional_effects(fit, resolution = 3L)),
                  c("xg", "z"))
  cep <- conditional_effects(fit, "xg", method = "predict", ndraws = 50L,
                             resolution = 4L)[[1L]]
  expect_true(all(is.finite(cep$estimate__)))
  # a level holds one value: within it there is no curve and no derivative
  nd <- data.frame(xg = c(-1, 0, 1), sg = 0.3, z = 0, g = factor("1"))
  expect_error(frm_curve(fit, newdata = nd, simultaneous = FALSE),
               "unique for each group")
  expect_error(frm_curve_deriv(fit, "xg", newdata = nd, simultaneous = FALSE),
               "no derivative or feature within a level")
  expect_error(frm_curve_feature(fit, "xg", type = "crossing", newdata = nd),
               "no derivative or feature within a level")
  # another covariate moves within a level, as before
  dz <- frm_curve_deriv(fit, "z", simultaneous = FALSE,
                        newdata = data.frame(xg = 0.5, sg = 0.3, z = c(0, 1),
                                             g = factor("1")))
  expect_equal(dz$.estimate, rep(unname(fixef(fit)["z", "Estimate"]), 2L),
               tolerance = 1e-6)
})

test_that("simulate(newdata) of me() beside cov = FALSE reads no response", {
  set.seed(8)
  d <- data.frame(t = 1:80, x = stats::rnorm(80), sx = 0.3)
  d$y <- 1 + 0.5 * d$x + as.numeric(stats::arima.sim(list(ar = 0.5), 80))
  fit <- frm(bf(y ~ me(x, sx) + ar(time = t, cov = FALSE)), data = d,
             family = gaussian())
  a <- data.frame(x = 0, sx = 0.3, t = 81:84, y = c(1, 2, 3, 4))
  b <- transform(a, y = c(9, -9, 9, -9))
  # the term starts from an empty past, so the response is not read; through
  # frm_linpred() the means moved with it, and the term counted twice
  expect_identical(simulate(fit, nsim = 20L, seed = 1, newdata = a),
                   simulate(fit, nsim = 20L, seed = 1, newdata = b))
  # the standard error of a cov = FALSE row holds this row's draw and, through
  # the residual, the last row's
  se <- frm_linpred(fit, newdata = a, se.fit = TRUE)$se.fit
  f0 <- frm_linpred(me_coef_only(fit), newdata = a, se.fit = TRUE)$se.fit
  expect_true(all(se > f0))
  fi <- fitted(fit, newdata = a)
  expect_equal(unname(fi[, "Est.Error"]), unname(se))
})

test_that("two grids at one x and one sdx read one draw in a difference", {
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  nd <- data.frame(x = c(-1, 1), sx = 0.5, z = 1)
  ct <- transform(nd, z = 0)
  cv <- frm_curve(fit, newdata = nd, contrast = ct, simultaneous = FALSE)
  # the coefficient of z, as a contrast in z reads it at one position
  expect_equal(cv$.se, rep(unname(fixef(fit)["z", "Est.Error"]), 2L),
               tolerance = 1e-6)
  # another x is another draw
  ct2 <- transform(nd, x = x + 1, z = 0)
  cv2 <- frm_curve(fit, newdata = nd, contrast = ct2, simultaneous = FALSE)
  expect_gt(cv2$.se[1L], 5 * cv$.se[1L])
  # the emmeans route takes a draw per row, as brms's draws do
  skip_if_not_installed("emmeans")
  em <- as.data.frame(pairs(emmeans::emmeans(fit, ~ z, at = list(
    z = c(0, 1), x = 1, sx = 0.5))))
  b1 <- fixef(fit)[grepl("^me", rownames(fixef(fit))), "Estimate"]
  expect_gt(em$SE, 0.95 * sqrt(2) * abs(b1) * 0.5)
})

test_that("marginaleffects warns of what its SE omits, and still answers", {
  skip_if_not_installed("marginaleffects")
  d <- mn_count()
  nd <- mn_grid_c()
  fit <- frm(bf(cnt ~ me(xo, sx) + z), family = poisson(), data = d)
  nz <- transform(nd, z = 0)
  expect_warning(marginaleffects::predictions(fit, newdata = nz),
                 class = "frmtmb_interop_se_part")
  expect_warning(marginaleffects::predictions(fit, newdata = nz), "se.fit")
  # a test of a coefficient is the right one, so a refusal stopped a right answer
  h <- allow_warnings(marginaleffects::hypotheses(fit, "b2 = 0"),
                      c(mfx_se_part, "positional indices"), mfx_se_part)
  cf <- fixef(fit)
  expect_equal(as.numeric(h$estimate), cf["z", "Estimate"])
  expect_equal(as.numeric(h$std.error), cf["z", "Est.Error"], tolerance = 1e-6)
  nd <- nz
  p <- marginaleffects::predictions(fit, newdata = nd, vcov = FALSE)
  # the estimate is the mean, as fitted() reads it
  expect_equal(p$estimate, unname(fitted(fit, newdata = nd)[, "Estimate"]))
  # the absent case: a fit with no such term gives no warning
  f0 <- frm(bf(cnt ~ xo), family = poisson(), data = d)
  expect_true(all(is.finite(
    marginaleffects::predictions(f0, newdata = nd)$std.error)))
})

test_that("the same warning reaches a field and a group-level term", {
  skip_if_not_installed("marginaleffects")
  set.seed(7411)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), g = factor(sample(1:12, n, TRUE)))
  d$y <- 1 + sin(2 * d$x) + stats::rnorm(12, 0, 0.8)[d$g] +
    stats::rnorm(n, 0, 0.5)
  nd <- data.frame(x = c(-1, 1), g = factor(c(1, 2), 1:12))
  # each of these holds a variance in no coefficient of get_coef()
  fs <- list(frm(y ~ s(x), data = d), frm(y ~ gp(x, k = 15), data = d))
  for (f in fs) {
    expect_warning(marginaleffects::predictions(f, newdata = nd),
                   class = "frmtmb_interop_se_part")
  }
  fg <- allow_warnings(frm(y ~ x + (1 | g), data = d), "boundary")
  p <- allow_warnings(marginaleffects::predictions(fg, newdata = nd),
                      mfx_se_part, "group-level effects")
  # what it warns of is true: the level's SE is not fitted()'s
  ft <- fitted(fg, newdata = nd)
  expect_gt(min(abs(p$std.error / ft[, "Est.Error"] - 1)), 0.1)
  # and with the group dropped it agrees, with no warning
  pn <- allow_warnings(
    marginaleffects::predictions(fg, newdata = nd, re_formula = NA),
    "not known to be supported")
  expect_equal(pn$std.error,
               unname(fitted(fg, newdata = nd, re_formula = NA)[, "Est.Error"]),
               tolerance = 1e-5)
})

test_that("fitted()'s interval holds when two predictors read the term", {
  # the draw's share read off the first predictor's gradient alone: me() in mu
  # and zi read 0.49 to 2.88 where the draw's quantiles are 1.30 to 2.02
  set.seed(8401)
  n <- 500
  d <- data.frame(lat = stats::rnorm(n), sx = 0.4)
  d$xo <- d$lat + stats::rnorm(n, 0, 0.4)
  eta <- 0.4 + 0.7 * d$lat
  d$cnt <- stats::rpois(n, exp(eta))
  d$zc <- ifelse(stats::runif(n) < stats::plogis(-1 + 0.8 * d$lat), 0L,
                 stats::rpois(n, 3))
  comp <- stats::rbinom(n, 1, 0.5)
  d$mixc <- stats::rpois(n, ifelse(comp == 1, exp(1.5 + 0.4 * d$lat),
                                   exp(0.2 + 0.3 * d$lat)))
  fits <- list(frm(bf(cnt ~ me(xo, sx)), family = poisson(), data = d),
               frm(bf(zc ~ me(xo, sx), zi ~ me(xo, sx)),
                   family = zero_inflated_poisson(), data = d),
               frm(bf(mixc ~ me(xo, sx)), family = mixture(poisson(), poisson()),
                   data = d))
  for (fit in fits) {
    for (sdx in c(0.4, 1.2)) {
      fv <- fitted(fit, newdata = data.frame(xo = 1, sx = sdx))
      # the expected response at the draw's own values, which frm_linpred() reads
      at <- data.frame(xo = 1 + sdx * stats::qnorm(stats::ppoints(2000)),
                       sx = sdx)
      v <- as.numeric(frm_linpred(me_coef_only(fit), newdata = at,
                                  type = "response"))
      q <- stats::quantile(v, c(0.025, 0.975), names = FALSE)
      # a zi that rises with the term turned this mean the other way
      expect_lt(fv[, "Q2.5"], fv[, "Q97.5"])
      # coefficient uncertainty widens the band a little past the draw's own
      expect_lt(abs(fv[, "Q2.5"] / q[1L] - 1), 0.1)
      expect_lt(abs(fv[, "Q97.5"] / q[2L] - 1), 0.1)
    }
  }
})

test_that("a category's interval stays inside 0 and 1, near the draw's", {
  d <- mn_count()
  d$ord <- as.integer(cut(0.4 + 0.7 * d$lat + stats::rlogis(nrow(d)),
                          c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
  fit <- frm(bf(ord ~ me(xo, sx)), family = cumulative(), data = d)
  fv <- fitted(fit, newdata = data.frame(xo = 1, sx = 1.2))
  at <- data.frame(xo = 1 + 1.2 * stats::qnorm(stats::ppoints(2000)), sx = 1.2)
  v <- frm_linpred(me_coef_only(fit), newdata = at, type = "response")
  # the two end categories are monotone in the draw
  for (k in c(1L, 4L)) {
    q <- stats::quantile(v[, k], c(0.025, 0.975), names = FALSE)
    expect_lt(abs(fv[1L, "Q2.5", k] / q[1L] - 1), 0.25)
    expect_lt(abs(fv[1L, "Q97.5", k] / q[2L] - 1), 0.25)
  }
  # a plain fit, far enough out for the Wald interval to pass 0
  f0 <- frm(bf(ord ~ xo), family = cumulative(), data = d)
  p0 <- fitted(f0, newdata = data.frame(xo = c(-6, 0, 6)))
  expect_true(all(p0[, "Q2.5", ] > 0 & p0[, "Q97.5", ] < 1))
  expect_true(all(p0[, "Q2.5", ] < p0[, "Estimate", ] &
                    p0[, "Estimate", ] < p0[, "Q97.5", ]))
})

test_that("a contrast reads one draw to a tolerance, not to the last bit", {
  # a grid built by arithmetic jumped by 20 times: 0.0485 at dx 0, 0.9979 at 1e-12
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  se_at <- function(dx) {
    frm_curve(fit, newdata = data.frame(x = 1, sx = 0.6, z = 1),
              contrast = data.frame(x = 1 + dx, sx = 0.6, z = 0),
              simultaneous = FALSE)$.se
  }
  expect_equal(se_at(1e-12), se_at(0))
  # two values of x past the tolerance are two draws
  expect_gt(se_at(1e-3), 5 * se_at(0))
})

test_that("one me() call is one latent across responses and dpars", {
  set.seed(7701)
  n <- 160
  lat <- stats::rnorm(n)
  d <- data.frame(xo = lat + stats::rnorm(n, 0, 0.4), sx = 0.4)
  d$y1 <- 1.5 * lat + stats::rnorm(n, 0, 0.4)
  d$y2 <- 1.2 * lat + stats::rnorm(n, 0, 0.4)
  fit <- frm(mvbf(bf(y1 ~ me(xo, sx)), bf(y2 ~ me(xo, sx)), rescor = FALSE),
             data = d, family = gaussian())
  nd <- data.frame(xo = c(0, 1), sx = c(1, 2))
  set.seed(3)
  dr <- predict(fit, newdata = nd, summary = FALSE, ndraws = 4000L,
                propagate_error = FALSE)
  s1 <- sigma(fit)[["y1"]]
  s2 <- sigma(fit)[["y2"]]
  b <- fixef(fit)[grepl("mexo", rownames(fixef(fit))), "Estimate"]
  target <- b[[1L]] * b[[2L]] * nd$sx^2 /
    sqrt(((b[[1L]] * nd$sx)^2 + s1^2) * ((b[[2L]] * nd$sx)^2 + s2^2))
  got <- vapply(1:2, function(i) stats::cor(dr[, i, 1L], dr[, i, 2L]), 0)
  expect_lt(max(abs(got - target)), 0.04)
  # one store: the two responses moved by one z
  sh <- fit
  sh[["me_draw"]] <- new.env(parent = emptyenv())
  set.seed(4)
  e1 <- frm_linpred(sh, newdata = nd, resp = "y1")
  e2 <- frm_linpred(sh, newdata = nd, resp = "y2")
  o1 <- frm_linpred(fit, newdata = nd, resp = "y1")
  o2 <- frm_linpred(fit, newdata = nd, resp = "y2")
  expect_equal(unname((e1 - o1) / b[[1L]]), unname((e2 - o2) / b[[2L]]))
})

test_that("predict(newdata) answers under cov = TRUE, as simulate() does", {
  set.seed(9)
  d <- data.frame(t = 1:120, x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + as.numeric(stats::arima.sim(list(ar = 0.5), 120))
  fit <- frm(bf(y ~ x + ar(time = t, cov = TRUE)), data = d,
             family = gaussian())
  # the term counts lags in the fitted time levels, as simulate() reads it
  nd <- data.frame(x = c(-1, 0, 1), t = 5:7)
  set.seed(2)
  p <- predict(fit, newdata = nd, summary = FALSE, ndraws = 3000L,
               propagate_error = FALSE)
  s <- simulate(fit, nsim = 3000L, seed = 2, newdata = nd)
  sm <- t(vapply(s, function(v) as.numeric(v), numeric(3L)))
  expect_lt(max(abs(colMeans(p) - colMeans(sm)) /
                  (apply(sm, 2, stats::sd) / sqrt(1500))), 4.5)
  expect_lt(max(abs(apply(p, 2, stats::sd) / apply(sm, 2, stats::sd) - 1)), 0.1)
  # the rows of one series correlate
  expect_gt(stats::cor(p[, 1L], p[, 2L]), 0.3)
})

test_that("a multivariate fit with no resp stacks, or refuses, and never picks", {
  set.seed(5)
  d <- data.frame(x = stats::rnorm(60))
  d$y <- d$x + stats::rnorm(60)
  d$y2 <- -d$x + stats::rnorm(60)
  fit <- frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE), data = d,
             family = gaussian())
  nd <- data.frame(x = c(-1, 1))
  m <- frm_linpred(fit, newdata = nd)
  expect_identical(colnames(m), c("y", "y2"))
  expect_equal(unname(m[, "y2"]),
               as.numeric(frm_linpred(fit, newdata = nd, resp = "y2")))
  p <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
  expect_identical(dim(p$se.fit), c(2L, 2L))
  expect_identical(dim(frm_linpred(fit, type = "response")), c(60L, 2L))
  expect_error(frm_lp_basis(fit, newdata = nd), "a design belongs to one")
  expect_error(frm_curve(fit, newdata = nd, simultaneous = FALSE),
               "a curve is read off one")
  expect_error(frm_curve_deriv(fit, "x", newdata = nd, simultaneous = FALSE),
               "a curve is read off one")
})

test_that("a row with no x or sdx reads NA with one warning, the others answer", {
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  whole <- fitted(fit, newdata = mn_grid())
  for (nd in list(transform(mn_grid(), sx = c(0.2, NA, 0.8)),
                  transform(mn_grid(), x = c(-1, NA, 1)))) {
    n <- 0L
    fv <- withCallingHandlers(fitted(fit, newdata = nd), warning = function(w) {
      if (inherits(w, "frmtmb_me_missing")) n <<- n + 1L
      invokeRestart("muffleWarning")
    })
    # one reader's many designs say it once
    expect_identical(n, 1L)
    expect_true(all(is.na(fv[2L, ])))
    expect_equal(fv[c(1L, 3L), "Estimate"], whole[c(1L, 3L), "Estimate"])
    lk <- allow_warnings(frm_linpred(fit, newdata = nd), "read NA", "read NA")
    expect_identical(unname(which(is.na(lk))), 2L)
    cv <- allow_warnings(frm_curve(fit, newdata = nd, simultaneous = FALSE),
                         "read NA", "read NA")
    expect_identical(which(is.na(cv$.estimate)), 2L)
    set.seed(5)
    sm <- allow_warnings(simulate(fit, newdata = nd, nsim = 2),
                         c("read NA", "not finite at the estimates"), "read NA")
    expect_identical(which(is.na(sm[, 1L])), 2L)
    expect_false(anyNA(sm[-2L, ]))
  }
  # a column that is absent is still an error, as no row has a value
  nd <- mn_grid()
  nd$sx <- NULL
  expect_error(fitted(fit, newdata = nd), "measurement SD")
})

test_that("fitted(newdata =) answers on the fit's own data with rows dropped", {
  d <- mn_data()
  d$x[c(5, 9)] <- NA
  fit <- suppressMessages(frm(bf(y ~ me(x, sx) + z), data = d,
                              family = gaussian()))
  fv <- allow_warnings(fitted(fit, newdata = d), "read NA", "read NA")
  expect_identical(nrow(fv), nrow(d))
  expect_identical(which(is.na(fv[, "Estimate"])), c(5L, 9L))
})

test_that("the curve says its band is not coefficient uncertainty alone", {
  fit <- frm(bf(y ~ me(x, sx) + z), data = mn_data(), family = gaussian())
  cv <- frm_curve(fit, newdata = mn_grid()[c("x", "sx", "z")],
                  simultaneous = FALSE)
  expect_output(print(cv), "not coefficient uncertainty")
  f0 <- frm(bf(y ~ x + z), data = mn_data(), family = gaussian())
  c0 <- frm_curve(f0, newdata = mn_grid()[c("x", "z")], simultaneous = FALSE)
  out <- utils::capture.output(print(c0))
  expect_false(any(grepl("not coefficient uncertainty", out)))
})

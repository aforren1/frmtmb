# mo() of an integer is coded from its minimum, as brms 2.23.0's
# get_mo_values() codes it: 1..4 has D = 3 and codes 0..3. Up to 0.75.0
# frmtmb coded the integers themselves, so 1..4 had D = 4 and a simplex
# cell for the category 0 that no row is in. Its maximum-likelihood
# coefficients were not identified (NaN standard errors on 3 of 4
# predictions) and frm_sample() sampled another model than brms. A
# negative integer was refused. Measured on 0.75.0 with seed 5:
# dev/mocode-repro.R, dev/mocode-findings.md.

mo_code_data <- function() {
  set.seed(5)
  n <- 300
  d <- data.frame(k = sample(1:4, n, replace = TRUE))
  d$k0 <- d$k - 1L
  d$k2 <- d$k + 1L
  d$kn <- d$k - 3L
  d$z <- rnorm(n)
  d$g <- factor(sample(1:8, n, replace = TRUE))
  d$y <- 2 + 4.5 * c(0, cumsum(c(0.6, 0.3, 0.1)))[d$k] + rnorm(n)
  d
}

mo_code_frame <- function(f, d) {
  fr <- frm(f, data = d, dry_run = "frame")
  fr$linpreds[[1L]]
}

test_that("an integer mo() is coded from its minimum", {
  d <- mo_code_data()
  for (v in c("k0", "k", "k2", "kn")) {
    lp <- mo_code_frame(stats::as.formula(paste0("y ~ mo(", v, ")")), d)
    mi <- lp$mo[[1L]]
    expect_identical(as.integer(mi$codes), d$k0, label = v)
    expect_identical(mi$D, 3L, label = v)
    expect_identical(mi$min, min(d[[v]]), label = v)
  }
  # a value between the minimum and the maximum that no row has is a
  # category all the same, as brms's D = max - min makes it
  dg <- data.frame(x = rep(c(1L, 2L, 5L), 10), y = rnorm(30))
  mi <- mo_code_frame(y ~ mo(x), dg)$mo[[1L]]
  expect_identical(mi$D, 4L)
  expect_identical(as.integer(mi$codes), dg$x - 1L)
})

test_that("mo() refuses what brms refuses and takes one or two values", {
  d <- data.frame(x = rep(c(1, 1.5, 2, 3), 5), y = rnorm(20))
  expect_error(frm(y ~ mo(x), data = d, dry_run = "frame"),
               "ordered factor or integers")
  # brms 2.23.0 builds two values as Jmo = 1 (a simplex of one weight)
  # and one value as Jmo = 0; through 0.76.0 frmtmb refused both ("mo()
  # needs at least 3 ordered categories")
  d$x <- rep(c(4L, 5L), 10)
  mi <- mo_code_frame(y ~ mo(x), d)$mo[[1L]]
  expect_identical(mi$D, 1L)
  expect_identical(as.integer(mi$codes), d$x - 4L)
  d$x <- 4L
  mi <- mo_code_frame(y ~ mo(x), d)$mo[[1L]]
  expect_identical(mi$D, 0L)
  expect_identical(as.integer(mi$codes), rep(0L, 20))
  d$x <- rep(c(4L, 5L, 6L, 6L), 5)
  expect_identical(mo_code_frame(y ~ mo(x), d)$mo[[1L]]$D, 2L)
})

mo_two_data <- function() {
  set.seed(5)
  n <- 200
  d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n))
  d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(n)
  d$x01 <- d$x - 1L
  d$xo <- factor(d$x, levels = 1:2, ordered = TRUE)
  d$x1 <- 4L
  d
}

test_that("mo() on two values is its coefficient times the code", {
  # brms's simplex[1] is fixed at 1, so mo(x) on two values is the linear
  # term b * (x - min(x)): the fit is lm()'s on the code, with the ML sigma
  d <- mo_two_data()
  fit <- frm(y ~ mo(x) + z, data = d)
  m <- stats::lm(y ~ x01 + z, data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(m)))
  # nlminb stops a few 1e-6 short of lm()'s exact solution here, as it
  # does on three categories, so the yardsticks are measured: the
  # estimates' own standard errors, and the gap between the ML and the
  # unbiased residual scale
  fe <- fixef(fit)
  se_lm <- sqrt(diag(stats::vcov(m)))[c(1, 3, 2)]
  s <- sqrt(m$df.residual / nrow(d))
  expect_lt(max(abs(fe[, "Estimate"] - coef(m)[c(1, 3, 2)]) / se_lm), 1e-3)
  expect_lt(max(abs(fe[, "Est.Error"] - se_lm * s)),
            1e-3 * max(abs(se_lm - se_lm * s)))
  # brms's summary row of the weight: 1, with no error
  sm <- summary(fit)$mo
  expect_identical(rownames(sm), "mox1[1]")
  expect_identical(unname(unlist(sm[1L, ])), c(1, 0, 1, 1))
  # an ordered factor with two levels is the same model
  fo <- frm(y ~ mo(xo) + z, data = d)
  expect_identical(as.numeric(logLik(fo)), as.numeric(logLik(fit)))
  # new data in the fitted range only, as brms's validate_newdata()
  nd <- data.frame(x = 1:2, z = 0)
  ex <- stats::predict(m, newdata = data.frame(x01 = 0:1, z = 0),
                       se.fit = TRUE)
  ft <- fitted(fit, newdata = nd)
  expect_lt(max(abs(ft[, "Estimate"] - ex$fit) / ex$se.fit), 1e-3)
  expect_lt(max(abs(ft[, "Est.Error"] - ex$se.fit * s)),
            1e-3 * max(abs(ex$se.fit - ex$se.fit * s)))
  expect_error(fitted(fit, newdata = data.frame(x = 3, z = 0)),
               "outside the fitted range 1..2")
  # brms's simo row has one weight, and its concentration has length 1
  expect_error(frm(y ~ mo(x) + z, data = d, prior = set_prior(
    "dirichlet(c(1, 1))", class = "simo", coef = "mox1")),
    "Expected input of length 1")
  fp <- frm(y ~ mo(x) + z, data = d, prior = set_prior(
    "dirichlet(2)", class = "simo", coef = "mox1"))
  expect_identical(as.numeric(logLik(fp)), as.numeric(logLik(fit)))
  # mo(x) * z takes it too
  fi <- frm(y ~ mo(x) * z, data = d)
  expect_equal(as.numeric(logLik(fi)),
               as.numeric(logLik(stats::lm(y ~ x01 * z, data = d))))
})

test_that("mo() on one value is 0 and leaves its coefficient open", {
  # brms 2.23.0 builds Jmo = 0; the term is 0 on every row, so the fit is
  # y ~ z, and the data do not determine the coefficient
  d <- mo_two_data()
  fit <- allow_warnings(frm(y ~ mo(x1) + z, data = d),
                        "Standard errors are not available",
                        require = "Standard errors are not available")
  m <- stats::lm(y ~ z, data = d)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(m)))
  expect_true(is.nan(fixef(fit)["mox1", "Est.Error"]))
  expect_null(summary(fit)$mo)
  nd <- data.frame(x1 = 4L, z = c(0, 1))
  ft <- fitted(fit, newdata = nd)
  expect_true(all(is.finite(ft[, "Est.Error"])))
  expect_error(fitted(fit, newdata = data.frame(x1 = 3L, z = 0)),
               "outside the fitted range 4..4")
  skip_if_not_installed("marginaleffects")
  # no prediction reads the coefficient, so marginaleffects' standard
  # errors are those of y ~ z; a NaN covariance made every one of them
  # NaN
  se <- marginaleffects::predictions(fit, newdata = nd)$std.error
  f0 <- frm(y ~ z, data = d)
  expect_true(all(is.finite(se)))
  expect_equal(se, marginaleffects::predictions(f0, newdata = nd)$std.error)
})

test_that("a one-value mo() term is no degree of freedom and no variance", {
  d <- mo_two_data()
  d$g <- factor(rep(1:8, length.out = nrow(d)))
  fit <- allow_warnings(frm(y ~ mo(x1) + z, data = d),
                        "Standard errors are not available")
  m <- stats::lm(y ~ z, data = d)
  # the term cannot change the likelihood, so it counts as lm() counts
  # an aliased column; the first build of this change gave df 4
  expect_equal(attr(logLik(fit), "df"), attr(logLik(m), "df"))
  expect_equal(AIC(fit), AIC(m))
  f0 <- frm(y ~ z, data = d)
  a <- anova(f0, fit)
  expect_identical(as.numeric(a[2L, "Chi Df"]), 0)
  # predict() draws the parameters: the lost coefficient is held at its
  # estimate, which no row reads, where its loss sent the draws to the
  # estimates alone, with a warning
  expect_no_warning(predict(fit, newdata = data.frame(x1 = 4L, z = 0),
                            ndraws = 50))
  # a group-level variation of a column that is 0 on every row is not a
  # parameter; fitted, it printed its start value 1 as an estimate
  expect_error(frm(y ~ mo(x1) + z + (mo(x1) | g), data = d),
               "x1 has one value")
})

test_that("brms codes one and two values as frmtmb does", {
  skip_unless_brms()
  d <- mo_two_data()
  for (v in c("x", "xo", "x1")) {
    f <- stats::as.formula(paste0("y ~ mo(", v, ") + z"))
    sd <- brms_standata(f, data = d, family = stats::gaussian())
    mi <- mo_code_frame(f, d)$mo[[1L]]
    expect_identical(as.integer(mi$codes), as.integer(sd$Xmo_1), label = v)
    expect_identical(mi$D, as.integer(sd$Jmo[1L]), label = v)
  }
})

test_that("mo() codes and simplex sizes are brms's Xmo and Jmo", {
  skip_unless_brms()
  d <- mo_code_data()
  d$kg <- c(1L, 2L, 5L)[(d$k %% 3L) + 1L]
  for (v in c("k0", "k", "k2", "kn", "kg")) {
    f <- stats::as.formula(paste0("y ~ mo(", v, ")"))
    sd <- brms_standata(f, data = d, family = stats::gaussian())
    mi <- mo_code_frame(f, d)$mo[[1L]]
    expect_identical(as.integer(mi$codes), as.integer(sd$Xmo_1), label = v)
    expect_identical(mi$D, as.integer(sd$Jmo[1L]), label = v)
  }
  # brms drops an ordered factor's unused levels before it codes them
  d$o <- factor(c("b", "c", "d")[(d$k %% 3L) + 1L],
                levels = c("a", "b", "c", "d"), ordered = TRUE)
  sd <- brms_standata(y ~ mo(o), data = d, family = stats::gaussian())
  mi <- mo_code_frame(y ~ mo(o), d)$mo[[1L]]
  expect_identical(as.integer(mi$codes), as.integer(sd$Xmo_1))
  expect_identical(mi$D, as.integer(sd$Jmo[1L]))
})

test_that("a shift of a mo() integer leaves the fit unchanged", {
  d <- mo_code_data()
  # an identity, not a tolerance: the four variables have the same codes,
  # so the four fits run the same tape from the same start
  for (rhs in c("mo(%s)", "mo(%s) * z", "mo(%s) + (mo(%s) | g)")) {
    fits <- lapply(c("k0", "k", "k2", "kn"), function(v) {
      frm(stats::as.formula(paste("y ~", gsub("%s", v, rhs, fixed = TRUE))),
          data = d)
    })
    ll <- vapply(fits, function(f) as.numeric(logLik(f)), 0)
    df <- vapply(fits, function(f) attr(logLik(f), "df"), 0)
    cf <- lapply(fits, function(f) unname(fixef(f)[, "Estimate"]))
    se <- lapply(fits, function(f) unname(fixef(f)[, "Est.Error"]))
    for (j in 2:4) {
      expect_identical(ll[j], ll[1L], label = rhs)
      expect_identical(df[j], df[1L], label = rhs)
      expect_identical(cf[[j]], cf[[1L]], label = rhs)
    }
    expect_true(all(is.finite(unlist(se))), label = rhs)
  }
})

test_that("new data are coded from the fitted minimum", {
  d <- mo_code_data()
  f1 <- frm(y ~ mo(k), data = d)
  f0 <- frm(y ~ mo(k0), data = d)
  p1 <- fitted(f1, newdata = data.frame(k = 1:4))
  p0 <- fitted(f0, newdata = data.frame(k0 = 0:3))
  expect_identical(p1, p0)
  expect_true(all(is.finite(p1[, "Est.Error"])))
  # a subset of the categories keeps its codes: brms codes new data
  # from the fitted data's minimum, not from their own
  expect_identical(fitted(f1, newdata = data.frame(k = 2:4)), p1[2:4, ])
  # brms refuses a value outside the fitted range (validate_newdata())
  expect_error(fitted(f1, newdata = data.frame(k = 0L)),
               "outside the fitted range 1..4")
  expect_error(fitted(f1, newdata = data.frame(k = 5L)),
               "outside the fitted range 1..4")
  fn <- frm(y ~ mo(kn), data = d)
  expect_identical(fitted(fn, newdata = data.frame(kn = -2:1)), p0)
  expect_error(predict(fn, newdata = data.frame(kn = -3L)),
               "outside the fitted range -2..1")
  # brms refuses a value that is not a whole number before the range:
  # 0.6, 2.5 and 4.4 were rounded to 1, 3 and 4 through 0.75.0
  for (v in c(0.6, 2.5, 4.4)) {
    expect_error(fitted(f1, newdata = data.frame(k = v)),
                 "defined at whole numbers only", label = v)
  }
  # each refusal names its own cause, and no base R warning escapes
  expect_error(fitted(f1, newdata = data.frame(k = NA_real_)),
               "new data have missing values")
  expect_error(fitted(f1, newdata = data.frame(
    k = factor(2, levels = 1:4, ordered = TRUE))),
    "fitted on integers, but new data give a factor")
  expect_error(fitted(f1, newdata = data.frame(k = "2")),
               "fitted on integers, but new data give a character vector")
  expect_error(frm(y ~ mo(mi(k)), data = d, dry_run = "frame"),
               "not the latent variable of mi[(][)] or me[(][)]")
})

test_that("a Dirichlet simo prior has the length brms gives it", {
  d <- mo_code_data()
  pr <- set_prior("dirichlet(c(2, 1, 1))", class = "simo", coef = "mok1")
  f1 <- frm(y ~ mo(k), data = d, prior = pr)
  pr0 <- set_prior("dirichlet(c(2, 1, 1))", class = "simo", coef = "mok01")
  f0 <- frm(y ~ mo(k0), data = d, prior = pr0)
  expect_identical(as.numeric(logLik(f1)), as.numeric(logLik(f0)))
  expect_error(frm(y ~ mo(k), data = d, prior = set_prior(
    "dirichlet(c(2, 1, 1, 1))", class = "simo", coef = "mok1")),
    "Expected input of length 3")
})

test_that("a frame coded from 0 by frmtmb 0.75.0 is refused when it differs", {
  d <- mo_code_data()
  f1 <- frm(y ~ mo(k), data = d)
  old <- f1
  # what 0.75.0 stored: the integers themselves, no minimum
  mi <- old$frame$linpreds[[1L]]$mo[[1L]]
  mi$min <- NULL
  mi$codes <- mi$codes + 1L
  old$frame$linpreds[[1L]]$mo[[1L]] <- mi
  expect_error(fitted(old), "frmtmb 0.75.0 or earlier")
  expect_error(fitted(old, newdata = data.frame(k = 2L)),
               "frmtmb 0.75.0 or earlier")
  # every method that reports an estimate, a df, a standard error or a
  # comparison refuses through one helper (mo_check_fit()). 0.75.0's
  # anova() of such a fit against its refit compared df 6 with df 5 at
  # one log-likelihood.
  for (m in c("fixef", "coef", "logLik", "AIC", "BIC", "vcov", "confint",
              "summary", "drop1")) {
    expect_error(do.call(m, list(old)), "frmtmb 0.75.0 or earlier",
                 label = m)
  }
  expect_error(anova(old, f1), "frmtmb 0.75.0 or earlier")
  # update() refits, and the refit is brms's model
  expect_identical(as.numeric(logLik(update(old))), as.numeric(logLik(f1)))
  # the absent case: on data that start at 0 the old coding is the same
  # model, and such a fit is kept
  f0 <- frm(y ~ mo(k0), data = d)
  old0 <- f0
  old0$frame$linpreds[[1L]]$mo[[1L]]$min <- NULL
  expect_identical(fitted(old0, newdata = data.frame(k0 = 0:3)),
                   fitted(f0, newdata = data.frame(k0 = 0:3)))
  expect_identical(fixef(old0), fixef(f0))
  expect_identical(logLik(old0), logLik(f0))
})

mo_gap_data <- function() {
  set.seed(7)
  d <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
  d$y <- d$x + rnorm(200)
  d
}

mo_gap_fit <- function() {
  allow_warnings(frm(y ~ mo(x), data = mo_gap_data()),
                 "Standard errors are not")
}

test_that("a prediction at a category no row is in has no standard error", {
  # x in 1, 2, 5: the categories 3 and 4 split two simplex steps that the
  # rows read only through their sum, so the data do not determine a
  # prediction there. Through 0.75.0 its standard error was finite
  # (0.0862 and 0.0951), because the delta method left the simplex out.
  fit <- mo_gap_fit()
  nd <- data.frame(x = 1:5)
  w <- "move along a direction the fit does not determine"
  p <- allow_warnings(fitted(fit, newdata = nd), w, require = w)
  expect_identical(is.nan(p[, "Est.Error"]), c(FALSE, FALSE, TRUE, TRUE,
                                               FALSE))
  s <- allow_warnings(frm_linpred(fit, newdata = nd, se.fit = TRUE), w,
                      require = w)
  expect_identical(unname(is.nan(s$se.fit)),
                   c(FALSE, FALSE, TRUE, TRUE, FALSE))
  ce <- allow_warnings(conditional_effects(fit, effects = "x"), w)[[1L]]
  expect_identical(ce$x, 1:5 + 0)
  expect_identical(is.nan(ce$se__), c(FALSE, FALSE, TRUE, TRUE, FALSE))
})

test_that("a prediction's standard error carries the simplex", {
  # the design of the delta method has d eta / d z for the simplex
  # coordinates z: compared with central differences of the prediction,
  # as a ratio to the largest entry
  d <- mo_code_data()
  cases <- list(
    list(f = y ~ mo(k), nd = data.frame(k = 1:4), re = NA),
    list(f = y ~ mo(k) * z, nd = data.frame(k = 1:4, z = c(-1, 0, 1, 2)),
         re = NA),
    list(f = y ~ mo(k) + (mo(k) | g),
         nd = data.frame(k = 1:4, g = factor(c(1, 3, 5, 8))), re = NULL))
  for (cs in cases) {
    fit <- frm(cs$f, data = d)
    lb <- frm_lp_basis(fit, newdata = cs$nd, re_formula = cs$re)
    zc <- which(get_joint_cov(fit)$names[lb$coef_pos] == "zeta1")
    expect_length(zc, 2L)
    z0 <- fit$estimates$zeta1
    h <- 1e-5
    num <- vapply(seq_along(z0), function(j) {
      at <- function(s) {
        f2 <- fit
        f2$estimates$zeta1[j] <- z0[j] + s * h
        frm_linpred(f2, newdata = cs$nd, re_formula = cs$re)
      }
      (at(1) - at(-1)) / (2 * h)
    }, numeric(nrow(cs$nd)))
    an <- as.matrix(lb$A)[, zc, drop = FALSE]
    expect_lt(max(abs(an - num)) / max(abs(num)), 1e-6,
              label = deparse1(cs$f))
  }
  # and on the fitted rows, where the group-level coefficient comes from
  # lp$Z
  fit <- frm(y ~ mo(k) + (mo(k) | g), data = d)
  lb <- frm_lp_basis(fit, re_formula = NULL)
  zc <- which(get_joint_cov(fit)$names[lb$coef_pos] == "zeta1")
  z0 <- fit$estimates$zeta1
  h <- 1e-5
  num <- vapply(seq_along(z0), function(j) {
    at <- function(s) {
      f2 <- fit
      f2$estimates$zeta1[j] <- z0[j] + s * h
      frm_linpred(f2, re_formula = NULL)
    }
    (at(1) - at(-1)) / (2 * h)
  }, numeric(nrow(d)))
  an <- as.matrix(lb$A)[, zc, drop = FALSE]
  expect_lt(max(abs(an - num)) / max(abs(num)), 1e-6)
})

test_that("conditional_effects() holds and moderates mo() as brms does", {
  d <- mo_code_data()
  fit <- frm(y ~ mo(k) + z, data = d)
  # brms 2.23.0 on these data holds k at round(median(k)) = 3 and
  # moderates with round(median +- mad) = 2, 3, 4 (dev/mocode-brms-fix1.R)
  ce <- conditional_effects(fit, effects = "z")[[1L]]
  expect_identical(unique(ce$k), round(stats::median(d$k)))
  p <- fitted(fit, newdata = data.frame(k = unique(ce$k), z = ce$z[1L]))
  expect_identical(ce$estimate__[1L], unname(p[1L, "Estimate"]))
  ce2 <- conditional_effects(fit, effects = "z:k")[[1L]]
  mk <- round(stats::median(d$k) + c(-1, 0, 1) * stats::mad(d$k))
  expect_identical(sort(unique(ce2$k)), mk)
  # the labels name the values the model was evaluated at
  expect_identical(sort(unique(as.numeric(as.character(ce2$effect2__)))),
                   mk)
})

test_that("emmeans refuses a mo() variable held between categories", {
  skip_if_not_installed("emmeans")
  d <- mo_code_data()
  fit <- frm(y ~ mo(k) + z, data = d)
  # brms 2.23.0 refuses too ("Invalid values in variable 'k'")
  expect_error(emmeans::emmeans(fit, ~ z),
               "emmeans holds 'k' at 2.486667, but mo[(]k[)] is defined")
  e <- emmeans::emmeans(fit, ~ z, at = list(k = 3))
  p <- frm_linpred(fit, newdata = data.frame(k = 3, z = mean(d$z)),
                   re_formula = NA)
  expect_equal(summary(e)$emmean, unname(p),
               tolerance = sqrt(.Machine$double.eps))
  # the remedy the refusal names rounds k alone: a bare function given
  # as cov.reduce would hold the continuous z at round(median(z)) = 0
  # instead of its mean
  expect_error(emmeans::emmeans(fit, ~ z),
               "cov.reduce = list[(]k = function[(]v[)] round")
  er <- emmeans::emmeans(fit, ~ z, cov.reduce = list(
    k = function(v) round(stats::median(v))))
  expect_identical(summary(er)$emmean, summary(e)$emmean)
})

test_that("marginaleffects' standard errors carry the mo() simplex", {
  skip_if_not_installed("marginaleffects")
  # with an interior simplex, y ~ mo(k) + z on four categories is
  # y ~ factor(k) + z in other coordinates, so the exact delta-method
  # standard error is lm()'s with the ML sigma. marginaleffects forms
  # J V J' over get_coef(), which through 0.75.0 left the simplex out:
  # 0.0846 and 0.1035 at k = 2 and 3, against 0.1220 and 0.1142
  d <- mo_code_data()
  fit <- frm(y ~ mo(k) + z, data = d)
  m <- stats::lm(y ~ factor(k) + z, data = d)
  nd <- data.frame(k = 1:4, z = 0)
  ex <- stats::predict(m, newdata = nd, se.fit = TRUE)$se.fit *
    sqrt(m$df.residual / nrow(d))
  se <- marginaleffects::predictions(fit, newdata = nd)$std.error
  # the standard error that leaves the simplex out, 0.75.0's, from the
  # basis and covariance of core's delta method: the yardstick for how
  # close marginaleffects must be
  lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
  jc <- get_joint_cov(fit)
  A <- as.matrix(lb$A)
  Vb <- as.matrix(jc$V)[lb$coef_pos, lb$coef_pos, drop = FALSE]
  zc <- which(jc$names[lb$coef_pos] == "zeta1")
  expect_length(zc, 2L)
  own <- frm_linpred(fit, newdata = nd, se.fit = TRUE)$se.fit
  expect_equal(unname(sqrt(rowSums((A %*% Vb) * A))), unname(own),
               tolerance = sqrt(.Machine$double.eps))
  A0 <- A[, -zc, drop = FALSE]
  se0 <- sqrt(rowSums((A0 %*% Vb[-zc, -zc, drop = FALSE]) * A0))
  # marginaleffects' Jacobian is a finite difference, about 1e-7
  # relative here; the simplex moves the standard error by 0.31
  expect_lt(max(abs(se / ex - 1)), 1e-4 * max(abs(se0 / ex - 1)))
  # a category no row is in: get_coef() carries the cumulative weights
  # at the observed categories, which the data determine, so the
  # standard error is NaN at the empty categories 3 and 4 alone, as
  # fitted() gives it. 0.75.0 gave finite ones there, and 0.76.0 NaN at
  # every category. At 1, 2 and 5 the exact one is lm()'s on factor(x)
  fg <- mo_gap_fit()
  dg <- mo_gap_data()
  mg <- stats::lm(y ~ factor(x), data = dg)
  exg <- stats::predict(mg, newdata = data.frame(x = c(1L, 2L, 5L)),
                        se.fit = TRUE)$se.fit *
    sqrt(mg$df.residual / nrow(dg))
  w <- "predictions read a mo() category that no row is in"
  pg <- allow_warnings(marginaleffects::predictions(
    fg, newdata = data.frame(x = 1:5)), w, require = w)
  expect_true(all(is.nan(pg$std.error[3:4])))
  # the yardstick as above, on the gap fit itself: how far the standard
  # error that leaves its simplex out is from the exact one
  ndo <- data.frame(x = c(1L, 2L, 5L))
  lbg <- allow_warnings(frm_lp_basis(fg, newdata = ndo, re_formula = NA),
                        "direction")
  jcg <- get_joint_cov(fg)
  Ag <- as.matrix(lbg$A)
  Vg <- as.matrix(jcg$V)[lbg$coef_pos, lbg$coef_pos, drop = FALSE]
  zg <- which(jcg$names[lbg$coef_pos] == "zeta1")
  Ag0 <- Ag[, -zg, drop = FALSE]
  se0g <- sqrt(rowSums((Ag0 %*% Vg[-zg, -zg, drop = FALSE]) * Ag0))
  yard <- max(abs(se0g / exg - 1))
  expect_gt(yard, 0)
  expect_lt(max(abs(pg$std.error[c(1, 2, 5)] / exg - 1)), 1e-4 * yard)
  expect_true(all(grepl("^zeta1_cum", setdiff(
    names(marginaleffects::get_coef(fg)),
    c("(Intercept)", "mox", "sigma_(Intercept)")))))
  # a contrast that reads no empty category is exact; one that does is
  # NaN
  ac <- marginaleffects::avg_comparisons(fg, variables = list(x = c(1, 5)))
  V <- stats::vcov(mg) * mg$df.residual / nrow(dg)
  expect_lt(abs(ac$std.error / sqrt(V[3, 3]) - 1), 1e-4 * yard)
  ag <- allow_warnings(marginaleffects::avg_comparisons(
    fg, variables = list(x = c(1, 3))), w)
  expect_true(is.nan(ag$std.error))
})

test_that("a gap fit that lost another direction warns in get_vcov()", {
  skip_if_not_installed("marginaleffects")
  # (mo(x) | g) with x in 1, 2, 5 stops at cor = 1 (seed 23 of the
  # review): the covariance on the determined directions is NaN itself.
  # The first build of this change returned NaN in every entry with no
  # warning; 0.76.0 kept the coefficient block and warned, and so does
  # this one
  set.seed(23)
  n <- 300
  d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                  g = factor(sample(1:10, n, TRUE)))
  d$y <- 1 + (0.4 + rnorm(10, 0, 0.2)[d$g]) * d$x + 0.5 * d$z + rnorm(n)
  fit <- allow_warnings(frm(y ~ mo(x) + z + (mo(x) | g), data = d),
                        "Optimizer did not report convergence")
  w <- "the fit lost another direction as well"
  V <- allow_warnings(marginaleffects::get_vcov(fit), w, require = w)
  cf <- c("(Intercept)", "z", "mox", "sigma_(Intercept)")
  expect_true(all(is.finite(V[cf, cf])))
  expect_true(all(is.nan(V["zeta1_cum1", ])))
})

test_that("dfbeta() keeps the coefficients' columns on a mo() fit", {
  # the simplex coordinates are in get_coef() for marginaleffects, and
  # not in the influence table, whose columns are the coefficients
  d <- mo_code_data()[1:40, ]
  fit <- frm(y ~ mo(k) + z, data = d)
  inf <- influence(fit)
  expect_identical(colnames(dfbeta(inf)),
                   c("(Intercept)", "z", "mok", "sigma_(Intercept)"))
  expect_identical(ncol(dfbetas(inf)), 4L)
  expect_length(cooks.distance(inf), 40L)
})

test_that("the grid of a mo() effect is the fitted term's categories", {
  # a subset() response codes its own rows: y's term was fitted on 2..4,
  # and its grid was every row's range 1..4 (0.75.0 stopped on a grid
  # of the wrong length)
  d <- mo_code_data()
  d$y2 <- d$y + stats::rnorm(nrow(d))
  d$s2 <- TRUE
  fit <- allow_warnings(
    frm(bf(y | subset(k > 1) ~ mo(k)) + bf(y2 | subset(s2) ~ mo(k)),
        data = d), "gradient")
  expect_identical(conditional_effects(fit, effects = "k",
                                       resp = "y")[[1L]]$k, c(2, 3, 4))
  expect_identical(conditional_effects(fit, effects = "k",
                                       resp = "y2")[[1L]]$k, c(1, 2, 3, 4))
})

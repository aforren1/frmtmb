# `.` in a formula, expanded against the data as brms 2.23.0 does
# (brms:::expand_dot_formula(): stats::terms(formula, data = data)), and
# a bf() delta in update(), read as brms:::update.brmsformula() reads it.

test_that("`.` expands to brms's formula and brms's design", {
  skip_on_cran()
  skip_if_not_installed("brms")
  dat <- data.frame(y = c(2, 5, 3, 8, 6, 9, 7, 12, 10, 11),
                    x1 = c(1, 3, 2, 5, 4, 6, 8, 7, 9, 10),
                    x2 = c(4, 1, 3, 2, 6, 5, 9, 7, 8, 10),
                    g = rep(1:2, 5))
  forms <- list(y ~ ., y ~ . - x2, y ~ 0 + ., y ~ . + (1 | g), y ~ x1 * .)
  for (f in forms) {
    fr <- frm(f, data = dat, dry_run = "frame")
    bt <- brms:::validate_formula(f, data = dat)$formula
    # the frame's location formula without its group-level terms
    expect_identical(deparse1(fr$spec$responses$y$dpars$mu$fixed),
                     deparse1(reformulas::RHSForm(reformulas::nobars(bt),
                                                  as.form = TRUE)),
                     label = deparse1(f))
    sd <- brms::standata(f, dat)
    cn <- colnames(fr$linpreds[["y.mu"]]$X)
    cn[cn == "(Intercept)"] <- "Intercept"
    expect_identical(cn, colnames(sd$X), label = deparse1(f))
  }
  # a parameter formula expands on its own, response included, as brms
  fr <- frm(bf(y ~ x1, sigma ~ .), data = dat, dry_run = "frame")
  sd <- brms::standata(brms::bf(y ~ x1, sigma ~ .), dat)
  cn <- colnames(fr$linpreds[["y.sigma"]]$X)
  cn[cn == "(Intercept)"] <- "Intercept"
  expect_identical(cn, colnames(sd$X_sigma))
  # an addition term's variable is not a predictor
  dw <- data.frame(y = dat$y, x1 = dat$x1, w = dat$x2 / 10)
  fr <- frm(y | weights(w) ~ ., data = dw, dry_run = "frame")
  expect_identical(colnames(fr$linpreds[["y.mu"]]$X), c("(Intercept)", "x1"))
})

test_that("a model with `.` is the model with the columns written out", {
  set.seed(3)
  d <- data.frame(y = rnorm(40), x1 = rnorm(40), x2 = rnorm(40))
  fit_dot <- frm(y ~ ., data = d)
  fit_out <- frm(y ~ x1 + x2, data = d)
  expect_identical(as.numeric(logLik(fit_dot)), as.numeric(logLik(fit_out)))
  expect_identical(fixef(fit_dot), fixef(fit_out))
  # the other entry points read the formula the same way
  expect_identical(get_prior(y ~ ., data = d), get_prior(y ~ x1 + x2, data = d))
  expect_identical(par_template(y ~ ., data = d),
                   par_template(y ~ x1 + x2, data = d))
  expect_identical(frm_simulate(y ~ ., d, nsim = 1, seed = 2,
                                newparams = par_template(fit_out)),
                   frm_simulate(y ~ x1 + x2, d, nsim = 1, seed = 2,
                                newparams = par_template(fit_out)))
})

test_that("update(newdata =) keeps the formula `.` expanded to", {
  # brms's update.brmsfit() reuses object$formula, which is the expanded
  # one, so new data with another column set fits the same model
  set.seed(6)
  d <- data.frame(y = rnorm(60), x1 = rnorm(60), x2 = rnorm(60))
  fit <- frm(y ~ ., data = d)
  # the call keeps the expanded formula, and prints as one
  expect_identical(deparse1(fit$call$formula), "y ~ x1 + x2")
  d2 <- d
  d2$extra <- rnorm(60)
  up <- update(fit, newdata = d2)
  expect_identical(rownames(fixef(up)), rownames(fixef(fit)))
  expect_identical(as.numeric(logLik(up)), as.numeric(logLik(fit)))
  expect_identical(as.numeric(logLik(update(fit, data = d2))),
                   as.numeric(logLik(fit)))
  # and new data without a column the model uses is refused, not fitted
  # without it
  expect_error(update(fit, newdata = d[, c("y", "x1")]), "x2",
               class = "frmtmb_error")
  # a bf() with `.` keeps its dpar formulas as written
  fb <- frm(bf(y ~ ., sigma ~ x1), data = d)
  upb <- update(fb, newdata = d2)
  expect_identical(as.numeric(logLik(upb)), as.numeric(logLik(fb)))
})

test_that("update() reads a bf() delta as brms does", {
  set.seed(4)
  d <- data.frame(x = rnorm(80), z = rnorm(80))
  d$y <- cut(d$x + rnorm(80), c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
  fit <- frm(y ~ x, data = d, family = cumulative())
  # brms's own spelling of a family change, update(fit3, bf(~ .,
  # family = acat())): the one-sided dot is the old formula
  up <- update(fit, bf(~ ., family = acat()))
  direct <- frm(y ~ x, data = d, family = acat())
  expect_identical(family(up)$family, "acat")
  expect_identical(as.numeric(logLik(up)), as.numeric(logLik(direct)))
  # a delta with a dpar formula keeps the location formula and adds it
  set.seed(5)
  dg <- data.frame(x = rnorm(80), z = rnorm(80))
  dg$y <- rnorm(80, 1 + dg$x, exp(0.3 * dg$z))
  fg <- frm(y ~ x, data = dg)
  up2 <- update(fg, bf(. ~ . + z, sigma ~ z))
  direct2 <- frm(bf(y ~ x + z, sigma ~ z), data = dg)
  expect_identical(as.numeric(logLik(up2)), as.numeric(logLik(direct2)))
})

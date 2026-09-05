## The rest of the frmtmb surface on a royston_parmar fit. A family is
## not finished when it fits; it is finished when the methods a user
## reaches for next either work or refuse for a stated reason.

sp_surface_fit <- function() {
  e <- new.env()
  utils::data("bc", package = "flexsurv", envir = e)
  bc <- e$bc
  bc$censored <- 1 - bc$censrec
  list(dat = bc,
       fit = frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                         family = royston_parmar(df = 2), data = bc))
}

test_that("summary, fixef and logLik work", {
  skip_if_not_installed("flexsurv")
  f <- sp_surface_fit()$fit
  expect_s3_class(f, "frmtmb_fit")
  expect_no_error(summary(f))
  fx <- frmtmb::fixef(f)
  expect_setequal(names(fx), c("mu", "gamma1", "gamma2"))
  expect_true(is.finite(as.numeric(stats::logLik(f))))
  expect_true(is.finite(stats::AIC(f)))
  expect_equal(stats::nobs(f), nrow(sp_surface_fit()$dat))
})

test_that("predict returns the spline coefficients, not a mean", {
  skip_if_not_installed("flexsurv")
  o <- sp_surface_fit()
  f <- o$fit
  p <- stats::predict(f, type = "link")
  expect_length(p, nrow(o$dat))
  expect_true(all(is.finite(p)))
  # every gamma is reachable by name
  for (dp in c("mu", "gamma1", "gamma2")) {
    v <- stats::predict(f, type = "link", dpar = dp)
    expect_true(all(is.finite(v)))
  }
  # the response scale has no meaning for this family and says so
  expect_error(stats::fitted(f), "royston_parmar")
  expect_error(stats::predict(f, type = "response"), "royston_parmar")
})

test_that("predict on newdata works and standard errors come with it", {
  skip_if_not_installed("flexsurv")
  f <- sp_surface_fit()$fit
  nd <- data.frame(group = factor(c("Good", "Poor"),
                                  levels = c("Good", "Medium", "Poor")))
  p <- stats::predict(f, newdata = nd, type = "link", se.fit = TRUE)
  expect_length(p$fit, 2L)
  expect_true(all(p$se.fit > 0))
  # a worse prognosis is a higher log cumulative hazard
  expect_gt(p$fit[2L], p$fit[1L])
})

test_that("residuals refuse for a stated reason", {
  skip_if_not_installed("flexsurv")
  f <- sp_surface_fit()$fit
  # "response" and "pearson" both need a fitted mean first and this
  # family declares none, so both refuse with the same message and
  # neither reaches the point of missing a variance function.
  # "deviance" is checked for a unit deviance before any mean is asked
  # for, so it refuses one step earlier and says so.
  expect_error(stats::residuals(f, type = "response"),
               "no mean on the response scale")
  expect_error(stats::residuals(f, type = "pearson"),
               "no mean on the response scale")
  expect_error(stats::residuals(f, type = "deviance"), "unit deviance")
})

test_that("par_template and set_prior reach the gamma dpars", {
  skip_if_not_installed("flexsurv")
  o <- sp_surface_fit()
  tmpl <- frmtmb::par_template(o$fit)
  nms <- names(unlist(tmpl))
  expect_true(any(grepl("gamma1_", nms, fixed = TRUE)))
  expect_true(any(grepl("gamma2_", nms, fixed = TRUE)))
  gp <- frmtmb::get_prior(frmtmb::bf(recyrs | cens(censored) ~ group),
                          family = royston_parmar(df = 2), data = o$dat)
  expect_true(is.data.frame(gp))
  expect_true(all(c("gamma1", "gamma2") %in% gp$dpar))
  # a prior on mu's slopes shrinks them
  f2 <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group),
                    family = royston_parmar(df = 2), data = o$dat,
                    prior = frmtmb::set_prior("normal(0, 0.05)", class = "b",
                                              dpar = "mu"))
  expect_lt(abs(unlist(frmtmb::fixef(f2))[["mu.groupPoor"]]),
            abs(unlist(frmtmb::fixef(o$fit))[["mu.groupPoor"]]))
})

test_that("a random effect and a smooth both reach a spline coefficient", {
  skip_if_not_installed("flexsurv")
  o <- sp_surface_fit()
  d <- o$dat
  set.seed(2)
  d$centre <- factor(rep(1:12, length.out = nrow(d)))
  fr <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group + (1 | centre)),
                    family = royston_parmar(df = 2), data = d)
  expect_s3_class(fr, "frmtmb_fit")
  expect_true(is.finite(as.numeric(stats::logLik(fr))))
  expect_gte(length(frmtmb::ranef(fr)), 1L)

  d$age <- stats::runif(nrow(d), 40, 80)
  fs <- frmtmb::frm(frmtmb::bf(recyrs | cens(censored) ~ group + s(age, k = 5)),
                    family = royston_parmar(df = 2), data = d)
  expect_s3_class(fs, "frmtmb_fit")
  # and the smooth it fitted is readable as a curve, with a band
  g <- data.frame(age = seq(45, 75, length.out = 20),
                  group = factor("Good", levels = levels(d$group)))
  cv <- frm_curve(fs, newdata = g, nsim = 2000, seed = 1)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-8)
  expect_true(all(is.finite(cv$.se)))
})

test_that("the compatibility rows this package registered resolve", {
  for (p in list(c("royston_parmar", "cens()"), c("royston_parmar", "trunc()"),
                 c("royston_parmar", "fitted"), c("frm_curve", "s()"),
                 c("frm_curve", "rr"), c("frm_curve", "predict"))) {
    r <- frmtmb::frm_compat(p[1L], p[2L])
    st <- if (is.data.frame(r)) r$status[1L] else r$status
    expect_true(st %in% c("works", "conditional", "refused", "untested",
                          "broken"),
                info = paste(p, collapse = " x "))
  }
  f <- frmtmb::frm_compat_features()
  expect_true(all(c("royston_parmar", "frm_curve") %in% f$key))
  expect_identical(frmtmb::frm_compat("royston_parmar", "fitted")$status,
                   "refused")
  # "conditional" and not "refused": re.form = NA works on an rr fit
  # and re.form = NULL does not, which test-curve.R measures
  expect_identical(frmtmb::frm_compat("frm_curve", "rr")$status,
                   "conditional")
})

test_that("print methods say what was checked", {
  o <- sp_surface_fit()
  set.seed(9)
  d <- data.frame(x = sort(stats::runif(150)))
  d$y <- sin(pi * d$x) + stats::rnorm(150, 0, 0.3)
  fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 6)),
                     family = stats::gaussian(), data = d)
  cv <- frm_curve(fit, newdata = data.frame(x = seq(0, 1, length.out = 10)),
                  nsim = 1000, seed = 1)
  expect_output(print(cv), "covariance checked against")
  expect_output(print(cv), "simultaneous")
  ft <- frm_curve_feature(fit, var = "x", type = "maximum",
                          newdata = data.frame(x = seq(0.05, 0.95,
                                                       length.out = 21)))
  expect_output(print(ft), "curve feature")
})

test_that("a smooth on gamma1 is a smooth time-varying effect", {
  skip_if_not_installed("flexsurv")
  o <- sp_surface_fit()
  d <- o$dat
  set.seed(3)
  d$age <- stats::runif(nrow(d), 40, 80)
  fs <- frmtmb::frm(
    frmtmb::bf(recyrs | cens(censored) ~ group, gamma1 ~ s(age, k = 5)),
    family = royston_parmar(df = 2), data = d)
  expect_s3_class(fs, "frmtmb_fit")
  expect_true(is.finite(as.numeric(stats::logLik(fs))))
  # and the fitted gamma1 surface is readable as a curve, which is the
  # case that first showed betad was missing from the design
  g <- data.frame(age = seq(45, 75, length.out = 15),
                  group = factor("Good", levels = levels(d$group)))
  cv <- frm_curve(fs, newdata = g, dpar = "gamma1", simultaneous = FALSE)
  expect_lt(attr(cv, "check")$cov_rel_error, 1e-10)
  expect_true(all(cv$.se > 0))
})

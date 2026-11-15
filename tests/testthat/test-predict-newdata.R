fit_np <- local({
  set.seed(31)
  n <- 400
  dd <- data.frame(
    x = rnorm(n),
    f = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    g = factor(rep(seq_len(20), length.out = n))
  )
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + c(a = 0, b = 1, c = -1)[dd$f] +
                  rnorm(20, 0, 0.8)[dd$g],
                exp(0.1 + 0.3 * dd$x))
  list(
    data = dd,
    fit = frm(bf(y ~ x + f + (1 | g), sigma ~ x) + gaussian(),
                 data = dd)
  )
})

test_that("newdata = training data reproduces in-sample predictions", {
  fit <- fit_np$fit; dd <- fit_np$data
  expect_equal(frm_linpred(fit, newdata = dd), frm_linpred(fit),
               tolerance = 1e-10)
  expect_equal(frm_linpred(fit, newdata = dd, type = "response"),
               frm_linpred(fit, type = "response"), tolerance = 1e-10)
  expect_equal(frm_linpred(fit, newdata = dd, dpar = "sigma",
                       type = "response"),
               frm_linpred(fit, dpar = "sigma", type = "response"),
               tolerance = 1e-10)
})

test_that("population-level predictions drop the random effects", {
  fit <- fit_np$fit; dd <- fit_np$data
  p0 <- frm_linpred(fit, newdata = dd, re_formula = NA)
  beta <- fixef_by_dpar(fit)$mu
  X <- model.matrix(~ x + f, dd)
  expect_equal(p0, unname(drop(X %*% beta)), tolerance = 1e-10,
               ignore_attr = TRUE)
})

test_that("factor levels and subsets round-trip", {
  fit <- fit_np$fit; dd <- fit_np$data
  nd <- dd[dd$f == "b", ][1:5, ]
  p_sub <- frm_linpred(fit, newdata = nd)
  expect_equal(p_sub, frm_linpred(fit)[which(dd$f == "b")[1:5]],
               tolerance = 1e-10)
})

test_that("new grouping levels error unless allowed", {
  fit <- fit_np$fit; dd <- fit_np$data
  nd <- dd[1:3, ]
  nd$g <- factor("999")
  expect_error(frm_linpred(fit, newdata = nd), "New levels")
  p_new <- frm_linpred(fit, newdata = nd, allow_new_levels = TRUE)
  p_pop <- frm_linpred(fit, newdata = nd, re_formula = NA)
  expect_equal(p_new, p_pop, tolerance = 1e-10)
})

test_that("se.fit matches glmmTMB delta-method standard errors", {
  skip_if_not_installed("glmmTMB")
  set.seed(32)
  n <- 300
  dd <- data.frame(x = rnorm(n),
                   g = factor(rep(seq_len(15), length.out = n)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(15, 0, 0.7)[dd$g], 1.2)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ref <- glmmTMB::glmmTMB(y ~ x + (1 | g), data = dd, REML = FALSE)

  nd <- data.frame(x = seq(-2, 2, length.out = 9),
                   g = factor(rep(1, 9), levels = levels(dd$g)))
  pf <- frm_linpred(fit, newdata = nd, re_formula = NA, se.fit = TRUE)
  # glmmTMB's own spelling, untouched: the rename is frmtmb's surface
  pr <- predict(ref, newdata = nd, re.form = NA, se.fit = TRUE)
  expect_vector_equal(pf$fit, pr$fit, tol = 1e-4)
  expect_vector_equal(pf$se.fit, pr$se.fit, tol = 1e-3)
})

test_that("response-scale se applies the chain rule", {
  set.seed(33)
  dd <- data.frame(x = rnorm(200))
  dd$y <- rpois(200, exp(0.5 + 0.4 * dd$x))
  fit <- frm(bf(y ~ x) + poisson(), data = dd)
  nd <- data.frame(x = c(-1, 0, 1))
  pl <- frm_linpred(fit, newdata = nd, se.fit = TRUE)
  pr <- frm_linpred(fit, newdata = nd, type = "response", se.fit = TRUE)
  expect_equal(pr$se.fit, exp(pl$fit) * pl$se.fit, tolerance = 1e-10)
})

test_that("predict(newdata =) draws a correlated residual on the new rows", {
  # the refusal said a draw follows the fitted rows, which is false of the route
  # simulate(newdata =) takes; test-me-newdata.R holds that route against it.
  # Only a family's own simulator still refuses (frmtmb.latent's
  # test-mixture-mvn.R)
  set.seed(3)
  da <- data.frame(g = factor(rep(1:10, each = 12)), t = rep(1:12, 10),
                   x = rnorm(120))
  da$y <- da$x + rnorm(120)
  fa <- frm(bf(y ~ x + ar(time = t, gr = g, cov = TRUE)),
            family = gaussian(), data = da)
  expect_identical(nrow(predict(fa, newdata = da[1:5, ])), 5L)
  # the guard absent: without the term, newdata is predicted
  fg <- frm(bf(y ~ x), family = gaussian(), data = da)
  expect_identical(nrow(predict(fg, newdata = da[1:5, ])), 5L)
})

test_that("a group-level design is rebuilt under the fit's contrasts", {
  # backlog 5353: the population-level design read the stored contrasts
  # and the group-level one the session's, so a fit made under contr.sum
  # refused every prediction once the option was reset ("columns:
  # (Intercept), fb, fc vs (Intercept), f1, f2")
  set.seed(7)
  dc <- data.frame(f = factor(sample(c("a", "b", "c"), 90, TRUE)),
                   g = factor(rep(1:10, each = 9)), x = stats::rnorm(90))
  dc$g2 <- factor(rep(c(2:10, 1), each = 9), levels = 1:10)
  dc$y <- as.numeric(dc$f) + stats::rnorm(10)[dc$g] + stats::rnorm(90)
  q <- function(...) suppressWarnings(suppressMessages(frm(...)))
  nd <- dc[c(1:6, 30), ]
  withr::with_options(list(contrasts = c("contr.sum", "contr.poly")), {
    fc <- q(y ~ f + (1 + f | g), data = dc)
    ref <- fitted(fc, newdata = nd)[, "Estimate"]
    ref_p <- fitted(fc, newdata = nd, re_formula = ~ (1 + f | g))[,
                                                                "Estimate"]
    fm <- q(y ~ x + (1 + f | mm(g, g2)), data = dc)
    ref_m <- fitted(fm, newdata = nd)[, "Estimate"]
  })
  expect_identical(fc$frame$re_blocks[[1L]]$cnms,
                   c("(Intercept)", "f1", "f2"))
  withr::with_options(list(contrasts = c("contr.treatment", "contr.poly")), {
    expect_equal(fitted(fc, newdata = nd)[, "Estimate"], ref)
    expect_equal(fitted(fc, newdata = nd,
                        re_formula = ~ (1 + f | g))[, "Estimate"], ref_p)
    expect_equal(fitted(fm, newdata = nd)[, "Estimate"], ref_m)
  })
  # the in-sample rows agree with newdata of the same rows
  expect_equal(unname(fitted(fc)[c(1:6, 30), "Estimate"]), unname(ref))
  # the guard's condition absent: a fit under the default contrasts
  # predicts the same under either option
  ft <- q(y ~ f + (1 + f | g), data = dc)
  r0 <- fitted(ft, newdata = nd)[, "Estimate"]
  withr::with_options(list(contrasts = c("contr.sum", "contr.poly")), {
    expect_equal(fitted(ft, newdata = nd)[, "Estimate"], r0)
  })
})

test_that("a transformed factor in a group-level term keeps the fit's contrasts", {
  # review of 2026-11-13, F1: `(1 + factor(cyl) | g)` stored no contrasts,
  # since the model frame holds `factor(cyl)` and no `cyl`, so the term was
  # rebuilt under the session's option: 0.483 off, silently, under
  # contr.helmert (whose columns are named as contr.sum's), and refused
  # under contr.treatment
  set.seed(7)
  dc <- data.frame(g = factor(rep(1:10, each = 9)),
                   cyl = sample(c(4, 6, 8), 90, TRUE))
  dc$y <- dc$cyl / 4 + stats::rnorm(10)[dc$g] + stats::rnorm(90)
  q <- function(...) suppressWarnings(suppressMessages(frm(...)))
  nd <- dc[c(1:6, 30), ]
  for (fml in list(y ~ 1 + (1 + factor(cyl) | g),
                   y ~ factor(cyl) + (1 + factor(cyl) | g))) {
    withr::with_options(list(contrasts = c("contr.sum", "contr.poly")), {
      fc <- q(fml, data = dc)
      ref <- fitted(fc, newdata = nd)[, "Estimate"]
    })
    expect_equal(unname(fitted(fc)[c(1:6, 30), "Estimate"]), unname(ref))
    for (o in list(c("contr.helmert", "contr.poly"),
                   c("contr.treatment", "contr.poly"))) {
      withr::with_options(list(contrasts = o), {
        expect_equal(fitted(fc, newdata = nd)[, "Estimate"], ref)
      })
    }
  }
})

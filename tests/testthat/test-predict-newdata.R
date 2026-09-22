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

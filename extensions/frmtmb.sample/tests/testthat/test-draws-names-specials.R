# The draws of a model with exp(x) carry brms's names. Through frmtmb
# 0.74.0 the parser's alias for exp() reached the design column, so the
# draws held b_.frm_expx, sd_g__.frm_expx and cor_g__Intercept__.frm_expx,
# and hypothesis(draws, "expx > 0") was refused.
#
# brms 2.23.0's variables() for this model, from a brmsfit built as
# help("brm") builds one from brms's own Stan code and standata
# (dev/brmsnames-1024-brms.R, case "re exp"; the b_, sd_ and cor_ names,
# in brms's order).

skip_on_cran()

test_that("draws of exp(x) are named as brms names them", {
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(1024)
  n <- 160
  d <- data.frame(x = stats::rnorm(n), g = factor(rep(1:16, 10)))
  u <- cbind(stats::rnorm(16, 0, 0.5), stats::rnorm(16, 0, 0.3))
  d$y <- 1 + (0.4 + u[d$g, 2]) * exp(d$x / 2) + u[d$g, 1] +
    stats::rnorm(n, 0, 0.5)
  fit <- frm(bf(y ~ exp(x) + (1 + exp(x) | g)), family = gaussian(),
             data = d)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 5)))
  v <- posterior::variables(ds)
  want <- c("b_Intercept", "b_expx", "sd_g__Intercept", "sd_g__expx",
            "cor_g__Intercept__expx")
  expect_identical(grep("^(b|sd|cor)_", v, value = TRUE), want)
  expect_false(any(grepl(".frm_", v, fixed = TRUE)))
  expect_identical(rownames(fixef(ds)), c("Intercept", "expx"))
  h <- hypothesis(ds, "expx > 0")
  expect_identical(nrow(h$hypothesis), 1L)
})

# A prior written with the design column's spelling or with brms's
# lands on the coefficient it names: under normal(3, 0.001) the draws of
# b_IzE2 sit at 3, and no other coefficient does. On frmtmb 0.74.0
# frm_sample() refused brms's spelling, "Prior target not found".
test_that("frm_sample() takes a coef in either spelling", {
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(31)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 1 + d$x + 0.3 * d$z^2 + stats::rnorm(n)
  fit <- frm(bf(y ~ x + I(z^2)), family = gaussian(), data = d)
  for (spell in c("I(z^2)", "IzE2")) {
    p <- set_prior("normal(3, 0.001)", class = "b", coef = spell)
    ds <- suppressWarnings(suppressMessages(
      frm_sample(fit, prior = p, chains = 1, iter = 200, refresh = 0,
                 seed = 5)))
    m <- colMeans(as.matrix(ds)[, c("b_x", "b_IzE2"), drop = FALSE])
    expect_identical(names(m)[abs(m - 3) < 0.01 * 3], "b_IzE2",
                     label = spell)
  }
})

# brms's fixef() of the draws is in brms's row order, which is the fit's
# fixef() order: every centered intercept first. brms 2.23.0 on this
# model gives Intercept, sigma_Intercept, x, sigma_expz
# (dev/brmsnames-1024-brms.R, case 32); the draws keep the template
# order, b_Intercept b_x b_sigma_Intercept b_sigma_expz, and fixef() of
# the draws followed it through frmtmb.sample 0.21.1.
test_that("fixef() of the draws is in brms's order", {
  skip_sampler()
  skip_if_not_installed("posterior")
  set.seed(32)
  n <- 160
  d <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(n, 0, exp(0.2 * exp(d$z / 2)))
  fit <- frm(bf(y ~ x, sigma ~ exp(z)), family = gaussian(), data = d)
  ds <- suppressWarnings(suppressMessages(
    frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 5)))
  want <- c("Intercept", "sigma_Intercept", "x", "sigma_expz")
  expect_identical(rownames(fixef(fit)), want)
  expect_identical(rownames(fixef(ds)), want)
  expect_identical(colnames(fixef(ds, summary = FALSE)), want)
})

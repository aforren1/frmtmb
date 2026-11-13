# integrated_design_zero(): which rows of new data have every design
# entry that multiplies an integrated value exactly 0. frmtmb.sample's
# laplace draws compute a prediction at such a row only; an entry that
# is not 0 refuses even where the entries cancel (review 2026-11-11 of
# lane sampfront, finding 1).

idz_data <- function(seed = 1530L) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(1:6, each = 8L)), t = rep(1:8, 6L),
                   x = stats::rnorm(48L))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(6L, 0, 0.6)[dd$g] +
    stats::rnorm(48L)
  dd
}

test_that("a row reads nothing integrated only where every entry is 0", {
  dd <- idz_data()
  q <- function(...) suppressWarnings(suppressMessages(frm(...)))
  nd <- data.frame(x = c(0, -1, 1, 0), g = factor(c(1, 2, 3, 9),
                                                  levels = 1:9))
  fs <- q(bf(y ~ x + (0 + x | g)), data = dd)
  expect_identical(integrated_design_zero(fs, nd),
                   c(TRUE, FALSE, FALSE, TRUE))
  # (1 + x | g) at x = -1 has entries 1 and -1, which cancel only when
  # the two effects are equal: not 0, so not skipped
  fi <- q(bf(y ~ x + (1 + x | g)), data = dd)
  expect_identical(integrated_design_zero(fi, nd), rep(FALSE, 4L))
  # a group-level term in sigma counts too, as every linear predictor of
  # the response does
  fsig <- q(bf(y ~ x + (0 + x | g), sigma ~ 1 + (1 | g)), data = dd)
  expect_identical(integrated_design_zero(fsig, nd), rep(FALSE, 4L))
  # the guard's condition absent: nothing integrated anywhere
  f0 <- q(bf(y ~ x), data = dd)
  expect_identical(integrated_design_zero(f0, nd), rep(TRUE, 4L))
  # a missing covariate is not a 0
  ndm <- data.frame(x = c(0, NA), g = factor(c(1, 2), levels = 1:6))
  expect_identical(integrated_design_zero(fs, ndm), c(TRUE, FALSE))
})

test_that("where the design does not say, no row is skipped", {
  dd <- idz_data()
  q <- function(...) suppressWarnings(suppressMessages(frm(...)))
  nd <- data.frame(x = 0, t = 3L, g = factor(1, levels = 1:6))
  # an autocorrelation term reads the fitted rows' residuals
  fa <- q(bf(y ~ x + (0 + x | g) + ar(time = t, gr = g)), data = dd)
  expect_false(integrated_design_zero(fa, nd))
  # a nonlinear body reads its parameters' designs through the body
  fn <- q(bf(y ~ a + b * x, a ~ 1, b ~ 0 + (1 | g), nl = TRUE), data = dd,
          prior = c(set_prior("normal(0, 2)", nlpar = "a")))
  expect_false(integrated_design_zero(fn, nd))
  # a smooth's basis is not 0 at a row inside its range
  fsm <- q(bf(y ~ s(x, k = 5)), data = dd)
  expect_false(integrated_design_zero(fsm, nd))
  # a design that does not build: a variable newdata lacks
  fs <- q(bf(y ~ x + (0 + x | g)), data = dd)
  expect_false(integrated_design_zero(fs, data.frame(g = factor(1))))
  expect_identical(integrated_design_zero(fs, data.frame()), logical(0))
})

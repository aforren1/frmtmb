# How the fit methods read `resp`, each against brms 2.23.0 as run in
# dev/predread-resp-brms-probe.R. On a
# univariate model brms ignores the argument whatever it holds. Seen to fail
# on 0.78.0 (rellib-r16): fitted(resp = character(0)) stopped with "attempt to
# select less than one element", resp = c("y", "y") with "Unknown response:
# 'yy'", and conditional_effects() of a multivariate fit returned the first
# response alone, or an internal "'length = 2' in coercion" for several.

rr_data <- function() {
  set.seed(5)
  n <- 80
  d <- data.frame(x = stats::rnorm(n))
  d$y <- d$x + stats::rnorm(n)
  d$y2 <- -d$x + stats::rnorm(n)
  d$y3 <- 0.5 * d$x + stats::rnorm(n)
  d
}

test_that("a univariate fit does not read `resp`", {
  d <- rr_data()
  fit <- frm(bf(y ~ x), data = d, family = gaussian())
  nd <- data.frame(x = c(-1, 1))
  want_f <- fitted(fit, newdata = nd)
  want_l <- frm_linpred(fit, newdata = nd)
  want_r <- residuals(fit)
  want_b <- frm_lp_basis(fit, newdata = nd)$eta
  want_c <- frm_curve(fit, newdata = nd, simultaneous = FALSE)$.estimate
  set.seed(1)
  want_p <- predict(fit, newdata = nd, ndraws = 20L)
  for (r in list(character(0), c("y", "y"), "y", "nosuch", c("a", "b"),
                 NA_character_, 1L)) {
    expect_identical(fitted(fit, newdata = nd, resp = r), want_f)
    expect_identical(frm_linpred(fit, newdata = nd, resp = r), want_l)
    expect_identical(residuals(fit, resp = r), want_r)
    expect_identical(frm_lp_basis(fit, newdata = nd, resp = r)$eta, want_b)
    expect_identical(frm_curve(fit, newdata = nd, resp = r,
                               simultaneous = FALSE)$.estimate, want_c)
    set.seed(1)
    expect_identical(predict(fit, newdata = nd, ndraws = 20L, resp = r),
                     want_p)
    expect_identical(names(conditional_effects(fit, "x", resp = r,
                                               resolution = 3L)), "x")
  }
})

test_that("a multivariate fit still reads `resp`", {
  d <- rr_data()
  fit <- frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE), data = d,
             family = gaussian())
  nd <- data.frame(x = c(-1, 1))
  expect_error(fitted(fit, newdata = nd, resp = "nosuch"), "Unknown response")
  expect_error(frm_linpred(fit, newdata = nd, resp = "nosuch"),
               "Unknown response")
  expect_error(predict(fit, newdata = nd, resp = "nosuch"), "Unknown response")
})

test_that("conditional_effects() gives one panel set per response", {
  d <- rr_data()
  fit <- frm(mvbf(bf(y ~ x), bf(y2 ~ x), bf(y3 ~ x), rescor = FALSE), data = d,
             family = gaussian())
  nm <- function(...) names(conditional_effects(fit, "x", ..., resolution = 3L))
  # names as brms 2.23.0 gives them, on one call each
  expect_identical(nm(), c("y.y_x", "y2.y2_x", "y3.y3_x"))
  expect_identical(nm(resp = character(0)), c("y.y_x", "y2.y2_x", "y3.y3_x"))
  expect_identical(nm(resp = "y2"), "y2.y2_x")
  expect_identical(nm(resp = c("y", "y3")), c("y.y_x", "y3.y3_x"))
  expect_identical(nm(resp = c("y3", "y")), c("y3.y3_x", "y.y_x"))
  expect_identical(nm(resp = c("y", "y3"), method = "predict", seed = 1),
                   c("y.y_x", "y3.y3_x"))
  expect_error(nm(resp = "nosuch"), "Invalid argument 'resp'")
  # brms's two panel sets of one response read as a typo here, the same
  # refusal draws give (draws_resp_set())
  expect_error(nm(resp = c("y", "y")), "twice")
  ce <- conditional_effects(fit, "x", resolution = 3L)
  expect_s3_class(ce, "frmtmb_conditional_effects")
  # each panel is its response's own curve
  for (r in c("y", "y2", "y3")) {
    p <- ce[[paste0(r, ".", r, "_x")]]
    expect_identical(attr(p, "response"), r)
    expect_equal(p$estimate__,
                 as.numeric(fitted(fit, newdata = data.frame(x = p$x),
                                   resp = r)[, "Estimate"]))
  }
  # print() and plot() read the joined object
  expect_output(print(ce[1L]), "x")
})

test_that("pp_mixture() takes one response of a multivariate mixture", {
  d <- rr_data()
  fm <- frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE), data = d,
            family = mixture(gaussian(), gaussian()),
            control = frmtmb_control(check_se = "ignore"))
  # brms's words for each (dev/predread-resp-brms-probe.R)
  for (r in list(NULL, character(0), c("y", "y"), c("y", "y2"))) {
    expect_error(pp_mixture(fm, resp = r), "must be a single variable name")
  }
  expect_error(pp_mixture(fm, resp = "nosuch"), "Invalid argument 'resp'")
  p1 <- pp_mixture(fm, resp = "y")
  p2 <- pp_mixture(fm, resp = "y2")
  expect_identical(dim(p1), c(nrow(d), 4L, 2L))
  expect_identical(dim(p2), dim(p1))
  expect_false(isTRUE(all.equal(p1, p2)))
  # a univariate mixture ignores it
  fu <- frm(bf(y ~ x), data = d, family = mixture(gaussian(), gaussian()),
            control = frmtmb_control(check_se = "ignore"))
  w <- pp_mixture(fu)
  expect_identical(pp_mixture(fu, resp = c("y", "y")), w)
  expect_identical(pp_mixture(fu, resp = character(0)), w)
  expect_identical(pp_mixture(fu, resp = "nosuch"), w)
})

test_that("the curve readers and emmeans do not read resp on one response", {
  d <- rr_data()
  fit <- frm(bf(y ~ x), data = d, family = gaussian())
  nd <- data.frame(x = seq(-2.1, 2, length.out = 9))
  dv <- frm_curve_deriv(fit, "x", newdata = nd, simultaneous = FALSE)$.estimate
  ft <- frm_curve_feature(fit, "x", type = "crossing",
                          newdata = nd)$.estimate
  em <- as.data.frame(emmeans::emmeans(fit, ~ x, at = list(x = 1)))$emmean
  gp <- marginaleffects::get_predict(fit, newdata = nd)$estimate
  for (r in list(character(0), c("y", "y"), "nosuch", NA_character_, 1L)) {
    expect_identical(frm_curve_deriv(fit, "x", newdata = nd, resp = r,
                                     simultaneous = FALSE)$.estimate, dv)
    expect_identical(frm_curve_feature(fit, "x", type = "crossing", newdata = nd,
                                       resp = r)$.estimate, ft)
    expect_identical(as.data.frame(emmeans::emmeans(
      fit, ~ x, at = list(x = 1), resp = r))$emmean, em)
    expect_identical(marginaleffects::get_predict(
      fit, newdata = nd, resp = r)$estimate, gp)
  }
})

test_that("panels are in brms's order: effects outside, responses inside", {
  d <- rr_data()
  d$f <- factor(rep(c("a", "b"), 40))
  fit <- frm(mvbf(bf(y2 ~ x), bf(y3 ~ x + f), rescor = FALSE), data = d,
             family = gaussian())
  # brms 2.23.0 gives y2.y2_f though y2 does not read f, a flat panel
  # (dev/review-predread-misc.R, block D)
  ce <- conditional_effects(fit, resolution = 3L)
  expect_named(ce, c("y2.y2_x", "y3.y3_x", "y2.y2_f", "y3.y3_f"))
  flat <- ce[["y2.y2_f"]]$estimate__
  expect_lt(diff(range(flat)), 1e-8 * max(1, abs(flat)))
  expect_named(conditional_effects(fit, c("x", "f"), resolution = 3L),
               c("y2.y2_x", "y3.y3_x", "y2.y2_f", "y3.y3_f"))
  expect_named(conditional_effects(fit, "x", resp = c("y3", "y2"),
                                   resolution = 3L),
               c("y3.y3_x", "y2.y2_x"))
})

test_that("a feature search and insight read a multivariate fit by response", {
  skip_if_not_installed("insight")
  set.seed(4)
  n <- 60
  d <- data.frame(y = stats::rnorm(n, 2), y2 = stats::rnorm(n, -1),
                  x = stats::rnorm(n))
  fit <- frm(bf(y ~ x) + bf(y2 ~ x) + set_rescor(FALSE), data = d,
             family = gaussian())
  nd <- data.frame(x = c(-1, 0, 1))
  # "'length = 2' in coercion to 'logical(1)'"
  expect_error(frm_curve_feature(fit, "x", newdata = nd), "Name it with resp =",
               class = "frmtmb_error")
  gp <- insight::get_predicted(fit, data = nd)
  expect_s3_class(gp, "get_predicted")
  expect_identical(names(gp), c("Row", "Response", "x", "Predicted"))
  expect_identical(gp$Response, rep(c("y", "y2"), each = 3L))
  expect_identical(gp$Row, rep(1:3, 2L))
  expect_equal(gp$Predicted, as.vector(frm_linpred(fit, newdata = nd)))
  # its interval is keyed to those two columns, which is how insight joins it
  full <- as.data.frame(gp)
  expect_identical(nrow(full), 6L)
  expect_true(all(c("SE", "CI_low", "CI_high") %in% names(full)))
})

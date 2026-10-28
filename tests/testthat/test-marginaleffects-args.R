# marginaleffects' own arguments on a frmtmb fit: `vcov =` reaches
# get_vcov(), and brms's `re_formula` and `allow_new_levels` reach
# get_predict(). Through 0.76.0 both methods dropped them, so a user's
# covariance, "HC3" or a cluster formula gave the model-based standard
# errors, and predictions(re_formula = NA) gave the conditional
# predictions, without a message. Measured on 0.76.0:
# dev/mointerop-me-probe.R (seed 5).

me_args_data <- function() {
  set.seed(5)
  n <- 200
  d <- data.frame(z = rnorm(n), g = factor(sample(1:8, n, TRUE)))
  d$y <- 1 + 0.5 * d$z + rnorm(8)[d$g] * 0.5 + rnorm(n)
  d$cnt <- rpois(n, exp(0.3 + 0.4 * d$z))
  d
}

test_that("vcov = takes a matrix, FALSE, a function and HC0 as for lm", {
  skip_if_not_installed("marginaleffects")
  d <- me_args_data()
  fit <- frm(y ~ z, data = d)
  m <- stats::lm(y ~ z, data = d)
  nd <- data.frame(z = c(-1, 0, 1))
  pr <- function(x, ...) {
    marginaleffects::predictions(x, newdata = nd, ...)$std.error
  }
  V <- marginaleffects::get_vcov(fit)
  # a doubled covariance doubles the variance, as it does for lm()
  r_frm <- pr(fit, vcov = 2 * V) / pr(fit)
  r_lm <- pr(m, vcov = 2 * stats::vcov(m)) / pr(m)
  expect_equal(r_frm, r_lm)
  expect_equal(pr(fit, vcov = function(x) 2 * marginaleffects::get_vcov(x)),
               pr(fit, vcov = 2 * V))
  expect_null(pr(fit, vcov = FALSE))
  expect_equal(pr(fit, vcov = TRUE), pr(fit))
  # in another order, by its names
  i <- rev(seq_len(nrow(V)))
  expect_equal(pr(fit, vcov = 2 * V[i, i]), pr(fit, vcov = 2 * V))
  # an unnamed one is read in get_coef()'s order, with marginaleffects'
  # warning
  w <- "no column names"
  expect_equal(allow_warnings(pr(fit, vcov = unname(2 * V)), w, require = w),
               pr(fit, vcov = 2 * V))
  expect_error(pr(fit, vcov = diag(2)), "get_coef[(][)] has 3 coefficients")
  # HC0 is the sandwich with each row its own cluster, which on a
  # gaussian linear model is White's HC0 exactly
  skip_if_not_installed("sandwich")
  expect_equal(pr(fit, vcov = "HC0"), pr(m, vcov = "HC0"))
})

test_that("vcov = a cluster formula is the glm cluster sandwich", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  d <- me_args_data()
  fit <- frm(cnt ~ z, data = d, family = stats::poisson())
  m <- stats::glm(cnt ~ z, data = d, family = stats::poisson())
  # frmtmb's estimates are nlminb's, a few 1e-7 from glm()'s, so the
  # yardstick is how far the sandwich is from the model-based covariance
  V0 <- unname(stats::vcov(m))
  for (v in list(~ g, "HC0")) {
    Vf <- unname(marginaleffects::get_vcov(fit, vcov = v))
    Vm <- unname(marginaleffects::get_vcov(m, vcov = v))
    expect_lt(max(abs(Vf - Vm)), 1e-3 * max(abs(Vm - V0)))
  }
  # the whole outer covariance of vcov_cluster(full = TRUE) is cut down
  # to get_coef()'s rows
  Vc <- vcov_cluster(fit, ~ g, "CR1", full = TRUE)
  expect_equal(marginaleffects::get_vcov(fit, vcov = Vc),
               marginaleffects::get_vcov(fit, vcov = ~ g))
})

test_that("vcov = strings frmtmb has no meaning for are refused by name", {
  skip_if_not_installed("marginaleffects")
  d <- me_args_data()
  fit <- frm(y ~ z, data = d)
  nd <- data.frame(z = 0)
  for (v in c("HC3", "HC2", "NeweyWest", "kenward-roger", "satterthwaite")) {
    expect_error(marginaleffects::predictions(fit, newdata = nd, vcov = v),
                 paste0("vcov = \"", v, "\" is not defined for a frmtmb fit"),
                 fixed = TRUE, label = v)
  }
  # a random effect that crosses the clusters has no per-row score
  fr <- frm(y ~ z + (1 | g), data = d)
  expect_error(marginaleffects::get_vcov(fr, vcov = "HC0"),
               "crosses `cluster`")
})

test_that("predictions(re_formula =) reaches the prediction", {
  skip_if_not_installed("marginaleffects")
  d <- me_args_data()
  fit <- frm(y ~ z + (1 | g), data = d)
  nd <- data.frame(z = c(-1, 0, 1), g = factor(1:3, levels = 1:8))
  # marginaleffects keeps a list of the arguments it knows per class,
  # with no way for a package to add its own, so it warns and passes the
  # argument on
  w <- "not known to be supported"
  p0 <- marginaleffects::predictions(fit, newdata = nd)
  pna <- allow_warnings(marginaleffects::predictions(fit, newdata = nd,
                                                     re_formula = NA),
                        w, require = w)
  p00 <- allow_warnings(marginaleffects::predictions(fit, newdata = nd,
                                                     re_formula = ~ 0),
                        w, require = w)
  pg <- allow_warnings(marginaleffects::predictions(fit, newdata = nd,
                                                    re_formula = ~ (1 | g)),
                       w, require = w)
  fna <- fitted(fit, newdata = nd, re_formula = NA)
  f0 <- fitted(fit, newdata = nd)
  expect_equal(pna$estimate, unname(fna[, "Estimate"]))
  expect_equal(p00$estimate, unname(fna[, "Estimate"]))
  expect_equal(p0$estimate, unname(f0[, "Estimate"]))
  expect_equal(pg$estimate, unname(f0[, "Estimate"]))
  # through 0.76.0 the argument was dropped: re_formula = NA returned the
  # conditional predictions
  expect_false(isTRUE(all.equal(pna$estimate, p0$estimate)))
  expect_true(all(is.finite(pna$std.error)))
})

test_that("predictions() keeps newdata's rows with unused factor levels", {
  skip_if_not_installed("marginaleffects")
  d <- me_args_data()
  fit <- frm(y ~ z + (1 | g), data = d)
  # marginaleffects pads new data with a row per factor level and drops
  # them by their negative rowid, which get_predict() renumbered: 3 rows
  # of a factor with 8 levels gave 11 predictions
  nd <- data.frame(z = c(-1, 0, 1), g = factor(1:3, levels = 1:8))
  p <- marginaleffects::predictions(fit, newdata = nd)
  expect_identical(nrow(p), 3L)
  expect_equal(p$estimate, unname(fitted(fit, newdata = nd)[, "Estimate"]))
})

# The sandwiches are formed on the directions the data determine
# (cluster_sandwich_full()). Formed from the covariance vcov() shows, a
# fit with an empty mo() category, a one-value mo() term or a variance
# at its boundary gave NaN for every standard error, silently. Measured:
# dev/mointerop-fix1.R.
me_gap_data <- function() {
  set.seed(7)
  d <- data.frame(x = sample(c(1L, 2L, 5L), 200, replace = TRUE))
  d$y <- d$x + rnorm(200)
  d$g <- factor(rep(1:20, each = 10))
  d
}

test_that("HC0, HC1 and ~ g are finite where the data determine them", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  d <- me_gap_data()
  fit <- allow_warnings(frm(y ~ mo(x), data = d), "Standard errors are not")
  m <- stats::lm(y ~ factor(x), data = d)
  nd <- data.frame(x = 1:5)
  ndo <- data.frame(x = c(1L, 2L, 5L))
  w <- "category that no row is in"
  pr <- function(v) {
    allow_warnings(marginaleffects::predictions(fit, newdata = nd,
                                                vcov = v)$std.error, w)
  }
  s0 <- pr(TRUE)
  for (v in c("HC0", "HC1")) {
    s <- pr(v)
    ref <- marginaleffects::predictions(m, newdata = ndo, vcov = v)$std.error
    expect_true(all(is.nan(s[3:4])), label = v)
    # the yardstick: how far the sandwich is from the model-based SE
    expect_lt(max(abs(s[c(1, 2, 5)] - ref)),
              1e-3 * max(abs(ref - s0[c(1, 2, 5)])), label = v)
  }
  sg <- pr(~ g)
  expect_true(all(is.finite(sg[c(1, 2, 5)])))
  expect_true(all(is.nan(sg[3:4])))
  # vcov_cluster() keeps NaN in the lost rows only
  Vc <- vcov_cluster(fit, ~ g, "CR1", full = TRUE)
  lost <- grepl("^zeta", rownames(Vc))
  expect_true(all(is.finite(Vc[!lost, !lost])))
  expect_true(all(is.nan(diag(Vc)[lost])))
  # a user's full matrix has NaN where the fit lost a parameter: warned
  wv <- "has no finite value for zeta1_cum1"
  allow_warnings(marginaleffects::get_vcov(
    fit, vcov = vcov(fit, full = TRUE)), wv, require = wv)
})

test_that("HC0 and HC1 on a one-value mo() fit are those of y ~ z", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  set.seed(5)
  n <- 200
  d <- data.frame(z = rnorm(n), x1 = 4L)
  d$y <- 1 + 0.5 * d$z + rnorm(n)
  fit <- allow_warnings(frm(y ~ mo(x1) + z, data = d),
                        "Standard errors are not")
  m <- stats::lm(y ~ z, data = d)
  nd <- data.frame(x1 = 4L, z = c(0, 1))
  s0 <- marginaleffects::predictions(fit, newdata = nd)$std.error
  for (v in c("HC0", "HC1")) {
    s <- marginaleffects::predictions(fit, newdata = nd, vcov = v)$std.error
    ref <- marginaleffects::predictions(m, newdata = nd[, "z", drop = FALSE],
                                        vcov = v)$std.error
    expect_lt(max(abs(s - ref)), 1e-3 * max(abs(ref - s0)), label = v)
  }
})

test_that("a cluster sandwich is finite on a fit with a variance at 0", {
  skip_if_not_installed("marginaleffects")
  set.seed(31)
  n <- 300
  d <- data.frame(x = rnorm(n), g = factor(sample(1:30, n, TRUE)),
                  w = runif(n, 0.5, 2), off = log(runif(n, 1, 3)))
  d$ypos <- rgamma(n, shape = 2, rate = 2 / exp(0.5 + 0.3 * d$x))
  d$ybin <- rbinom(n, 1, pnorm(0.2 + 0.6 * d$x))
  d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + d$off))
  d$y <- 1 + 0.5 * d$x + rnorm(n) * (1 + abs(d$x))
  fit <- allow_warnings(frm(y ~ x + (1 | g), data = d), "boundary")
  V <- marginaleffects::get_vcov(fit, vcov = ~ g)
  expect_true(all(is.finite(V)))
  # the sandwich, not the default: 0.76.0 ignored vcov =
  expect_false(isTRUE(all.equal(V, marginaleffects::get_vcov(fit))))
  # vcov_cluster() itself: NaN only in the lost row (every entry before)
  Vc <- vcov_cluster(fit, ~ g, "CR1")
  expect_true(all(is.finite(Vc)))
})

test_that("HC1 and stata are lm()'s; the other strings give their reason", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  d <- me_args_data()
  fit <- frm(y ~ z, data = d)
  m <- stats::lm(y ~ z, data = d)
  nd <- data.frame(z = c(-1, 0, 1))
  s0 <- marginaleffects::predictions(fit, newdata = nd)$std.error
  for (v in c("HC1", "stata")) {
    s <- marginaleffects::predictions(fit, newdata = nd, vcov = v)$std.error
    ref <- marginaleffects::predictions(m, newdata = nd, vcov = v)$std.error
    expect_lt(max(abs(s - ref)), 1e-3 * max(abs(ref - s0)), label = v)
  }
  expect_error(marginaleffects::get_vcov(fit, vcov = "HC3"), "leverage")
  expect_error(marginaleffects::get_vcov(fit, vcov = "robust"),
               "maps \"robust\" to \"HC3\"", fixed = TRUE)
  expect_error(marginaleffects::get_vcov(fit, vcov = "HAC"), "time order")
  fr <- frm(y ~ z + (1 | g), data = d, REML = TRUE)
  expect_error(marginaleffects::get_vcov(fr, vcov = vcov(fr, full = TRUE)),
               "this fit is REML")
})

test_that("get_predict() refuses by name the arguments it does not take", {
  skip_if_not_installed("marginaleffects")
  d <- me_args_data()
  fit <- frm(y ~ z + (1 | g), data = d)
  nd <- data.frame(z = c(-1, 0, 1), g = factor(1:3, levels = 1:8))
  # lme4's re.form = NA returned the conditional predictions
  args <- list(re.form = NA, allow.new.levels = TRUE,
               sample_new_levels = "gaussian", ndraws = 10, draw_ids = 1:3,
               zitype = "response")
  for (a in names(args)) {
    cl <- c(list(fit, newdata = nd), args[a])
    expect_error(allow_warnings(do.call(marginaleffects::predictions, cl),
                                "not known to be supported"),
                 paste0("passed ", a, " ="), fixed = TRUE, label = a)
  }
  # dpar = reaches the prediction
  fs <- frm(bf(y ~ z, sigma ~ z), data = d)
  ndz <- data.frame(z = c(-1, 0, 1))
  p <- allow_warnings(marginaleffects::predictions(fs, newdata = ndz,
                                                   dpar = "sigma"),
                      "not known to be supported")
  expect_equal(p$estimate,
               unname(fitted(fs, newdata = ndz, dpar = "sigma")[, 1]))
})

test_that("the gap warning counts the predictions asked for", {
  skip_if_not_installed("marginaleffects")
  set.seed(51)
  n <- 600
  d <- data.frame(x = sample(c(1, 2, 5), n, TRUE), z = rnorm(n),
                  g = factor(sample(1:20, n, TRUE)))
  d$y <- 1 + rnorm(20, 0, 0.6)[d$g] + 0.4 * (d$x - 1) + 0.5 * d$z +
    rnorm(n)
  fit <- allow_warnings(frm(y ~ mo(x) + z + (1 | g), data = d),
                        "Standard errors are not")
  # marginaleffects pads the 4 rows with one per level of g; the warning
  # said "1 of 24"
  w <- "1 of 4 predictions read a mo() category"
  allow_warnings(marginaleffects::predictions(
    fit, newdata = data.frame(x = c(1, 2, 3, 5), z = 0,
                              g = factor(1, levels = 1:20))), w, require = w)
})

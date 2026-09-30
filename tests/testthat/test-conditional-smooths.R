# conditional_smooths(): brms's display of each smooth term on its own.
# The values are checked against quantities the fit reports through
# its public surface (a contrast of frm_linpred(), the delta method over
# frm_lp_basis()), not against the engine's own internals; the layout
# against brms is in dev/postfit2-brms-compare.R.

cs_data <- function(n = 150, seed = 1) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n), z = stats::runif(n),
                  f = factor(sample(c("a", "b"), n, TRUE)),
                  g = factor(sample(letters[1:5], n, TRUE)))
  d$y <- sin(2 * pi * d$x) + d$z * (d$f == "b") + stats::rnorm(n, 0, 0.3)
  d
}

test_that("conditional_smooths() draws each term with brms's keys and grid", {
  d <- cs_data()
  fit <- frm(bf(y ~ f + s(x) + s(z, by = f)), family = gaussian(), data = d)
  cs <- conditional_smooths(fit)
  expect_s3_class(cs, "frmtmb_conditional_effects")
  expect_true(attr(cs, "smooths_only"))
  expect_identical(names(cs), c("mu: s(x)", "mu: s(z,by=f)"))
  # brms's grid: 100 points of x; z by the two levels of f
  expect_equal(nrow(cs[[1]]), 100)
  expect_equal(nrow(cs[[2]]), 200)
  expect_identical(names(cs[[1]]),
                   c("x", "effect1__", "cond__", "estimate__", "se__",
                     "lower__", "upper__"))
  expect_identical(names(cs[[2]]),
                   c("z", "f", "effect1__", "effect2__", "cond__",
                     "estimate__", "se__", "lower__", "upper__"))
  expect_identical(attr(cs[[2]], "effects"), c("z", "f"))
  expect_identical(attr(cs[[1]], "response"), "mu: s(x)")
  expect_false(attr(cs[[1]], "surface"))
  expect_equal(nrow(attr(cs[[1]], "points")), nrow(d))
  # whitespace in the selector does not matter, as in brms
  expect_identical(names(conditional_smooths(fit, smooths = "s(z, by = f)")),
                   "mu: s(z,by=f)")
})

test_that("the drawn term is the fitted term and its band the delta method", {
  d <- cs_data()
  fit <- frm(bf(y ~ z + s(x)), family = gaussian(), data = d)
  cs <- conditional_smooths(fit, resolution = 25)[[1]]
  # the term up to its constant: a contrast of the linear predictor
  # over a grid that varies x alone
  nd <- data.frame(x = cs$x, z = 0.3)
  eta <- frm_linpred(fit, newdata = nd, type = "link")
  # absolute gaps, bounded by a fraction of the curve's own scale
  expect_lt(max(abs((cs$estimate__ - cs$estimate__[1]) - (eta - eta[1]))),
            1e-8 * max(abs(eta)))
  # and its constant: mgcv centers the term over the data
  at_data <- conditional_smooths(fit, int_conditions = list(x = d$x))[[1]]
  expect_lt(abs(mean(at_data$estimate__)),
            1e-8 * max(abs(at_data$estimate__)))
  # the delta method with the intercept and z columns removed
  lb <- frm_lp_basis(fit, newdata = nd, re_formula = NA)
  keep <- !lb$coef_names %in% c("beta.(Intercept)", "beta.z")
  A <- lb$A[, keep, drop = FALSE]
  se <- sqrt(rowSums((A %*% lb$V[keep, keep]) * A))
  expect_lt(max(abs(cs$se__ - se)), 1e-8 * max(se))
  expect_lt(max(abs((cs$upper__ - cs$estimate__) -
                      stats::qnorm(0.975) * cs$se__)), 1e-8 * max(se))
})

test_that("conditional_smooths() refuses as brms does, and without draws", {
  d <- cs_data()
  fit <- frm(bf(y ~ s(x)), family = gaussian(), data = d)
  expect_error(conditional_smooths(fit, smooths = "s3"),
               "No valid smooth terms found in the model")
  lin <- frm(bf(y ~ x), family = gaussian(), data = d)
  expect_error(conditional_smooths(lin),
               "No valid smooth terms found in the model")
  for (a in list(list(spaghetti = TRUE), list(ndraws = 10),
                 list(draw_ids = 1:3))) {
    expect_error(do.call(conditional_smooths, c(list(fit), a)),
                 paste0("`", names(a), "` has nothing to act on"))
  }
  expect_error(conditional_smooths(fit, nonesuch = 1), "nonesuch")
})

test_that("a grouping-factor smooth is drawn, as re_formula = NA keeps it", {
  d <- cs_data()
  fit <- allow_warnings(
    frm(bf(y ~ s(x, g, bs = "fs", k = 5)), family = gaussian(), data = d),
    c("convergence", "Hessian", "singular"))
  cs <- conditional_smooths(fit)
  expect_identical(names(cs), "mu: s(x,g,bs=\"fs\",k=5)")
  expect_identical(levels(cs[[1]]$g), levels(d$g))
  expect_equal(nrow(cs[[1]]), 100 * nlevels(d$g))
  expect_true(all(is.finite(cs[[1]]$estimate__)))
})

test_that("a two-covariate term is a surface, and too_far trims it", {
  d <- cs_data()
  fit <- frm(bf(y ~ t2(x, z)), family = gaussian(), data = d)
  full <- conditional_smooths(fit, resolution = 12)[[1]]
  expect_true(attr(full, "surface"))
  expect_equal(nrow(full), 144)
  expect_true(is.numeric(full$effect2__))
  far <- mgcv::exclude.too.far(full$x, full$z, d$x, d$z, dist = 0.1)
  trimmed <- conditional_smooths(fit, resolution = 12, too_far = 0.1)[[1]]
  expect_equal(nrow(trimmed), sum(!far))
  expect_gt(sum(far), 0)
  expect_equal(trimmed$estimate__, full$estimate__[!far])
  flat <- conditional_smooths(fit, resolution = 12, surface = FALSE)[[1]]
  expect_false(attr(flat, "surface"))
  expect_equal(nrow(flat), 36)
  expect_s3_class(flat$effect2__, "factor")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(conditional_smooths(fit, resolution = 12,
                                         too_far = 0.1), ask = FALSE))
})

test_that("keys name the parameter, the response and the nonlinear parameter", {
  d <- cs_data()
  d$y2 <- d$y + stats::rnorm(nrow(d))
  fs <- frm(bf(y ~ s(x), sigma ~ s(z)), family = gaussian(), data = d)
  expect_identical(names(conditional_smooths(fs)),
                   c("mu: s(x)", "sigma: s(z)"))
  fm <- frm(mvbf(bf(y ~ s(x)), bf(y2 ~ s(z))), data = d)
  expect_identical(names(conditional_smooths(fm)),
                   c("mu_y: s(x)", "mu_y2: s(z)"))
  fn <- frm(bf(y ~ a + b * z, a ~ s(x), b ~ 1, nl = TRUE),
            family = gaussian(), data = d)
  expect_identical(names(conditional_smooths(fn)), "mu_a: s(x)")
})

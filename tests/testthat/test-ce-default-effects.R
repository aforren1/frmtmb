# conditional_effects()'s default displays are brms's get_all_effects():
# one display per term of every parameter's predictor, of the variables
# that term reads. The expected names are brms 2.23.0's, recorded on the
# same formulas (dev/parity2-brms-probe.R, dev/parity2-findings.md).

ce_shape_data <- function() {
  set.seed(1)
  n <- 150
  d <- data.frame(x = stats::runif(n, 0.1, 2), z = stats::rnorm(n),
                  k = sample(1:4, n, TRUE), w = stats::rnorm(n),
                  g = factor(sample(letters[1:8], n, TRUE)))
  d$y <- 1 + d$x * d$z + 0.3 * d$k + stats::rnorm(n)
  d$yo <- cut(d$y, stats::quantile(d$y, c(0, 0.3, 0.6, 1)),
              include.lowest = TRUE, labels = FALSE)
  d$y2 <- 0.5 * d$x + stats::rnorm(n)
  d$x2 <- d$x + stats::rnorm(n, sd = 0.3)
  d$x2[seq(3, n, 10)] <- NA
  d
}

ce_default_names <- function(f, family = gaussian(), data,
                             allowed = character()) {
  fit <- suppressMessages(allow_warnings(
    frm(f, family = family, data = data), allowed))
  names(suppressMessages(conditional_effects(fit, resolution = 5)))
}

test_that("default effects follow brms's get_all_effects() term by term", {
  d <- ce_shape_data()
  # brms 2.23.0: x:z, then x z x:z, k z k:z, x z, x w
  expect_identical(ce_default_names(bf(y ~ exp(x):z), data = d), "x:z")
  expect_identical(ce_default_names(bf(y ~ t2(x, z)), data = d),
                   c("x", "z", "x:z"))
  expect_identical(ce_default_names(bf(y ~ mo(k):z), data = d),
                   c("k", "z", "k:z"))
  # every distributional parameter's predictor, whichever is drawn
  expect_identical(ce_default_names(bf(y ~ x, sigma ~ exp(z)), data = d),
                   c("x", "z"))
  # a variable only a group-level term reads is not a display. Under
  # OpenBLAS (0.3.26 and 0.3.32 with LAPACK) this fit stops at nlminb's
  # code 8 at a measured maximum, which frmtmb reports; the names are the
  # point here
  expect_identical(ce_default_names(bf(y ~ x + (1 + z | g), sigma ~ w),
                                    data = d,
                                    allowed = "A Newton step from the stop"),
                   c("x", "w"))
  # a term of three variables is left out; its pairs come from x * z * w
  expect_identical(ce_default_names(bf(y ~ x * z * w), data = d),
                   c("x", "z", "w", "x:z", "x:w", "z:w"))
  expect_identical(ce_default_names(bf(y ~ poly(x, 2) + z + x:w), data = d),
                   c("x", "z", "x:w"))
  expect_identical(ce_default_names(bf(y ~ a * exp(b * x), a ~ 1 + w,
                                       b ~ z, nl = TRUE), data = d),
                   c("x", "w", "z"))
  # a pair puts its numeric variable first, as brms's
  # prepare_cond_data(reorder = TRUE) does: brms 2.23.0 gives f x x:f
  d$f <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  expect_identical(ce_default_names(bf(y ~ f * x), data = d),
                   c("f", "x", "x:f"))
  # effects the caller names keep their order, as in brms
  fit <- frm(bf(y ~ f * x), family = gaussian(), data = d)
  expect_identical(names(conditional_effects(fit, effects = "f:x",
                                             resolution = 5)), "f:x")
  # a numeric grouping variable counts as a factor only where re_formula
  # keeps its term, as brms reads its group variables off
  # update_re_terms(formula, re_formula): brms 2.23.0 gives kk x kk:x at
  # re_formula = NA and kk x x:kk at re_formula = NULL
  set.seed(21)
  dk <- data.frame(x = stats::rnorm(200), kk = sample(1:6, 200, TRUE))
  dk$y <- 1 + 0.5 * dk$x + 0.2 * dk$kk + stats::rnorm(6)[dk$kk] +
    stats::rnorm(200)
  fk <- frm(bf(y ~ kk * x + (1 | kk)), family = gaussian(), data = dk)
  expect_identical(names(conditional_effects(fk, resolution = 5)),
                   c("kk", "x", "kk:x"))
  expect_identical(names(conditional_effects(fk, resolution = 5,
                                             re_formula = NULL)),
                   c("kk", "x", "x:kk"))
})

test_that("default effects of a cs() term and of a multivariate model", {
  d <- ce_shape_data()
  # brms 2.23.0: z x (the fixed effect, then the cs() term); the display
  # here is per category, so each name carries brms's cats__ key
  fit <- suppressMessages(frm(bf(yo ~ cs(x) + z), family = sratio(),
                              data = d))
  ce <- suppressMessages(conditional_effects(fit, resolution = 5))
  expect_identical(names(ce), c("z:cats__", "x:cats__"))
  # brms 2.23.0: y.y_x y2.y2_x y.y_z y2.y2_z
  expect_identical(ce_default_names(bf(y ~ x) + bf(y2 ~ z) +
                                      set_rescor(FALSE), data = d),
                   c("y.y_x", "y2.y2_x", "y.y_z", "y2.y2_z"))
  # every response's displays are drawn for one response too, as brms
  # reads them off the whole model
  fit2 <- suppressMessages(frm(bf(y ~ x) + bf(y2 ~ z) + set_rescor(FALSE),
                               family = gaussian(), data = d))
  ce2 <- suppressMessages(conditional_effects(fit2, resp = "y2",
                                              resolution = 5))
  expect_identical(names(ce2), c("y2.y2_x", "y2.y2_z"))
})

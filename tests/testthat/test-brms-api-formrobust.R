# Small brms spellings frmtmb refused (dev/formrobust-findings.md
# section 6): the deprecated lower-case `0 + intercept` (ledger rows
# standata:970, :974), frm(drop_unused_levels = ) (standata:1133),
# autocorrelation terms as objects (brm:112) and autocor() on a fit
# (brmsfit-methods:112, :113).

ba_data <- local({
  set.seed(61)
  d <- data.frame(x = stats::rnorm(40), g = gl(4, 10))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(40)
  d
})

ba_warn <- "Reserved variable name 'intercept' is deprecated"

test_that("0 + intercept is a column of ones, with brms's warning", {
  fit <- NULL
  allow_warnings(fit <- frm(bf(y ~ 0 + intercept + x), data = ba_data),
                 ba_warn, require = ba_warn)
  # brms keeps the name: b_intercept
  expect_identical(rownames(fixef(fit)), c("intercept", "x"))
  ref <- frm(bf(y ~ x), data = ba_data)
  # the model y ~ 1 + x, reached by a different start
  expect_lt(abs(as.numeric(logLik(fit) - logLik(ref))),
            1e-10 * abs(as.numeric(logLik(ref))))
  expect_lt(abs(fixef(fit)["intercept", 1] - fixef(ref)["Intercept", 1]) /
              fixef(ref)["Intercept", 2], 1e-4)
  # newdata needs no intercept column, as brms fills it
  nd <- ba_data[1:3, "x", drop = FALSE]
  expect_identical(dim(fitted(fit, newdata = nd)), c(3L, 4L))
  # the guard-absent case: in a formula with an intercept, `intercept`
  # is an ordinary variable and nothing warns
  d <- ba_data
  d$intercept <- stats::rnorm(40)
  f2 <- frm(bf(y ~ intercept + x), data = d)
  expect_identical(rownames(fixef(f2)), c("Intercept", "intercept", "x"))
})

test_that("a data column intercept that is not all ones is refused", {
  d <- ba_data
  d$intercept <- 2
  allow_warnings(
    expect_error(frm(bf(y ~ 0 + intercept + x), data = d),
                 "Variable name 'intercept' is reserved", fixed = TRUE),
    ba_warn)
})

test_that("brms's standata has the same design", {
  skip_unless_brms()
  sd <- suppressWarnings(brms_standata(y ~ 0 + intercept + x,
                                       data = ba_data))
  fr <- NULL
  allow_warnings(fr <- frm(bf(y ~ 0 + intercept + x), data = ba_data,
                           dry_run = "frame"), ba_warn)
  X <- fr$linpreds[[1L]]$X
  expect_identical(colnames(X), colnames(sd$X))
  expect_identical(matrix(as.numeric(X), nrow(X)),
                   matrix(as.numeric(sd$X), nrow(sd$X)))
})

test_that("frm(drop_unused_levels = ) is brms's argument", {
  dd <- data.frame(y = c(0.2, -0.4, 1.1, 0.3, -0.8, 0.5, 0.9, -0.1, 0.4,
                         0.7),
                   x = factor(c("a", "b"), levels = c("a", "b", "c")))
  fr1 <- frm(bf(y ~ x), data = dd, dry_run = "frame")
  expect_false("c" %in% levels(fr1$data_frame$x))
  fr0 <- NULL
  expect_message(
    fr0 <- frm(bf(y ~ x), data = dd, drop_unused_levels = FALSE,
               dry_run = "frame"),
    "dropping column(s): xc", fixed = TRUE)
  expect_identical(levels(fr0$data_frame$x), c("a", "b", "c"))
  expect_error(frm(bf(y ~ x), data = dd, drop_unused_levels = NA),
               "drop_unused_levels", fixed = TRUE)
})

test_that("autocorrelation terms exist as objects, and join by acformula", {
  expect_s3_class(ar(week, subj), "frmtmb_ac_term")
  expect_identical(arma(t, g, p = 2)$label, "arma(t, g, p = 2)")
  expect_identical(frmtmb::cosy(t, g)$label, "cosy(t, g)")
  # brms refuses the object added to a formula, in these words
  expect_error(bf(y ~ 1) + arma(x),
               "Autocorrelation terms can only be specified", fixed = TRUE)
  expect_error(bf(y ~ 1) + bf(z ~ 1) + ma(x),
               "not added to a 'mvbrmsformula' object", fixed = TRUE)
  expect_identical(deparse1((bf(y ~ x) + acformula(~ ar(t, g)))$formula),
                   "y ~ x + ar(t, g)")
  expect_identical(deparse1(bf(y ~ x, autocor = ~ ar(t, g))$formula),
                   "y ~ x + ar(t, g)")
  expect_error(acformula(~ x), "must contain at least one autocorrelation",
               fixed = TRUE)
  expect_error(bf(y ~ 1) + bf(z ~ 1) + acformula(~ ar(t, g)),
               "without the response variable name", fixed = TRUE)
  set.seed(62)
  dg <- expand.grid(t = 1:5, g = factor(1:12))
  dg$y <- stats::rnorm(60)
  f0 <- frm(bf(y ~ 1 + ar(t, g)), data = dg)
  f1 <- frm(bf(y ~ 1) + acformula(~ ar(t, g)), data = dg)
  expect_identical(logLik(f0), logLik(f1))
})

test_that("autocor() is brms's deprecated accessor", {
  set.seed(63)
  dg <- expand.grid(t = 1:5, g = factor(1:12))
  dg$y <- stats::rnorm(60)
  fit <- frm(bf(y ~ 1 + ar(t, g)), data = dg)
  msg <- "Method 'autocor' is deprecated"
  out <- 1
  allow_warnings(out <- autocor(fit), msg, require = msg)
  expect_null(out)
  allow_warnings(expect_error(autocor(fit, resp = "z"), "Invalid argument"),
                 msg)
})

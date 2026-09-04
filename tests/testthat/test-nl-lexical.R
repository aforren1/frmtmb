# Which names of a bf(nl = TRUE) body are columns of `data`, and which
# are objects of the formula environment.
#
# Every symbol of a nonlinear body that is not a nonlinear parameter is
# a request for a column. `drop_nl_lexical_datavars()` (R/frame.R) takes
# back that request for a name that is not in `data` and that resolves
# in the formula environment to something `model.frame()` could never
# hold as a column: a function, a list (a data.frame is one), an
# environment, a language object. Those are the arguments of a helper -
# `solve_pk(pk_dyn, ..., events = doses)` - not column references.
#
# The boundary matters in both directions, so this file pins both
# sides: the refused types resolve lexically, and vectors and matrices
# do not. A vector in the environment IS a legal model-frame variable
# and `model.frame()` already finds it there; the tests below freeze
# that behavior rather than change it.

lex_dd <- function(n = 30, seed = 7) {
  set.seed(seed)
  data.frame(y = stats::rnorm(n), x = stats::rnorm(n),
             g = factor(rep_len(letters[1:5], n)))
}

# --- the refused types resolve lexically ----------------------------

test_that("an env data.frame passed to a helper is not a data column", {
  d <- lex_dd()
  tab <- data.frame(k = c(2, 3, 5))
  scale_by <- function(df, x) x * nrow(df)
  form <- bf(y ~ b0 * scale_by(tab, x), b0 ~ 1, nl = TRUE)

  fr <- frm(form + gaussian(), data = d, dry_run = "frame",
            start = list(beta = 1))
  expect_false("tab" %in% names(fr$linpreds[["y.mu"]]$data_list))
  expect_setequal(names(fr$linpreds[["y.mu"]]$data_list), "x")
  # `scale_by` is the head of the call, so all.vars() never asked for it
  expect_identical(fr$linpreds[["y.mu"]]$nl_lexical, "tab")

  # and the fit is the fit of the same body written inline
  f_sym <- frm(form + gaussian(), data = d, start = list(beta = 1))
  f_inl <- frm(bf(y ~ b0 * scale_by(data.frame(k = c(2, 3, 5)), x),
                  b0 ~ 1, nl = TRUE) + gaussian(),
               data = d, start = list(beta = 1))
  expect_equal(as.numeric(logLik(f_sym)), as.numeric(logLik(f_inl)),
               tolerance = 1e-12)
})

test_that("an env list, environment and formula are not data columns", {
  d <- lex_dd()
  cfg <- list(k = 2)                      # a plain list
  box <- new.env()                        # an environment
  assign("k", 3, envir = box)
  fo <- ~ 1                               # a language object
  use <- function(l, e, f, x) x * l$k * get("k", envir = e) * length(f)
  form <- bf(y ~ b0 * use(cfg, box, fo, x), b0 ~ 1, nl = TRUE)

  fr <- frm(form + gaussian(), data = d, dry_run = "frame",
            start = list(beta = 1))
  expect_setequal(names(fr$linpreds[["y.mu"]]$data_list), "x")
  expect_setequal(fr$linpreds[["y.mu"]]$nl_lexical,
                  c("cfg", "box", "fo"))

  f <- frm(form + gaussian(), data = d, start = list(beta = 1))
  expect_s3_class(f, "frmtmb_fit")
})

# --- a column still wins --------------------------------------------

test_that("a data column wins over a same-named environment object", {
  d <- lex_dd()
  d$tab <- d$x                            # a real column called `tab`
  tab <- data.frame(k = 1:3)              # ... and an object of that name
  form <- bf(y ~ b0 * tab, b0 ~ 1, nl = TRUE)

  fr <- frm(form + gaussian(), data = d, dry_run = "frame",
            start = list(beta = 1))
  expect_true("tab" %in% names(fr$linpreds[["y.mu"]]$data_list))
  expect_identical(fr$linpreds[["y.mu"]]$nl_lexical, character(0))
  expect_equal(fr$linpreds[["y.mu"]]$data_list$tab, d$x)
})

test_that("a list column named in a body still fails at the frame", {
  # the G2.12 contract: an unusable column is refused by name, and the
  # same-named data.frame in the environment does not rescue it
  d <- lex_dd()
  d$tab <- replicate(nrow(d), list(1:2), simplify = FALSE)
  tab <- data.frame(k = 1:3)
  form <- bf(y ~ b0 * tab, b0 ~ 1, nl = TRUE)
  expect_error(frm(form + gaussian(), data = d, dry_run = "frame",
                   start = list(beta = 1)),
               "invalid type \\(list\\).*tab")

  # and the unused list column is still ignored
  form2 <- bf(y ~ b0 * x, b0 ~ 1, nl = TRUE)
  expect_s3_class(frm(form2 + gaussian(), data = d, dry_run = "frame",
                      start = list(beta = 1)),
                  "frmtmb_frame")
})

# --- the other side of the boundary, frozen -------------------------

test_that("an env vector or matrix still reaches the model frame", {
  # model.frame() resolves a variable through the formula environment
  # when `data` has no such column, and a vector or a matrix is a legal
  # column. The exemption must not take those: a matrix covariate is a
  # feature, and a stray vector of the right length has to keep
  # behaving as it does for a linear formula.
  d <- lex_dd()
  z <- stats::rnorm(nrow(d))
  fr <- frm(bf(y ~ b0 * z, b0 ~ 1, nl = TRUE) + gaussian(), data = d,
            dry_run = "frame", start = list(beta = 1))
  expect_true("z" %in% names(fr$linpreds[["y.mu"]]$data_list))
  expect_identical(fr$linpreds[["y.mu"]]$nl_lexical, character(0))

  m <- matrix(stats::rnorm(2 * nrow(d)), nrow(d), 2)
  fr_m <- frm(bf(y ~ b0 * m[, 1], b0 ~ 1, nl = TRUE) + gaussian(),
              data = d, dry_run = "frame", start = list(beta = 1))
  expect_true("m" %in% names(fr_m$linpreds[["y.mu"]]$data_list))
  expect_true(is.matrix(fr_m$linpreds[["y.mu"]]$data_list$m))
})

test_that("a bad env vector and a missing name fail as they always did", {
  d <- lex_dd()
  short <- stats::rnorm(3)
  expect_error(frm(bf(y ~ b0 * short, b0 ~ 1, nl = TRUE) + gaussian(),
                   data = d, dry_run = "frame", start = list(beta = 1)),
               "variable lengths differ")
  expect_error(frm(bf(y ~ b0 * no_such_thing, b0 ~ 1, nl = TRUE) +
                     gaussian(),
                   data = d, dry_run = "frame", start = list(beta = 1)),
               "object 'no_such_thing' not found")
})

test_that("nl_lexical_only draws the line at model.frame's types", {
  # the refused types
  expect_true(frmtmb:::nl_lexical_only(sum))
  expect_true(frmtmb:::nl_lexical_only(data.frame(a = 1)))
  expect_true(frmtmb:::nl_lexical_only(list(1, 2)))
  expect_true(frmtmb:::nl_lexical_only(globalenv()))
  expect_true(frmtmb:::nl_lexical_only(~ x))
  expect_true(frmtmb:::nl_lexical_only(quote(a + b)))
  # the legal columns
  expect_false(frmtmb:::nl_lexical_only(1:3))
  expect_false(frmtmb:::nl_lexical_only(c(1.5, 2.5)))
  expect_false(frmtmb:::nl_lexical_only(letters))
  expect_false(frmtmb:::nl_lexical_only(factor("a")))
  expect_false(frmtmb:::nl_lexical_only(Sys.Date()))
  expect_false(frmtmb:::nl_lexical_only(matrix(1:4, 2)))
  # an absent name: model.frame() must keep reporting it
  expect_false(frmtmb:::nl_lexical_only(NULL))
})

# Tabular input forms beyond plain data.frame. The whole pipeline runs
# through stats::model.frame(), which normalizes these; the tests pin
# that the normalization is lossless.

#' @srrstats {G2.7} Tibble and data.table inputs are accepted and give
#'   bit-identical fits to the same data as a plain data.frame; tested
#'   below on a random-intercept model.
#' @srrstats {G2.10} Single columns are extracted internally with `[[`
#'   and through `stats::model.frame()`, never with drop-sensitive
#'   single-bracket indexing, so tabular classes that change the default
#'   drop behavior (tibbles) fit identically. The tibble test below is
#'   the evidence: a drop= divergence would change the design matrix.
#' @srrstats {G2.11} Columns carrying extra classes or attributes
#'   (units-style vectors) pass through `model.frame()` with their
#'   numeric content intact; tested below against the unclassed fit.
#' @srrstats {G2.12} List columns are handled: an unused list column is
#'   ignored, and a list column named in the formula fails at
#'   `model.frame()` with an error naming the variable and its invalid
#'   type, before any estimation runs.
NULL

tab_dd <- local({
  set.seed(41)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
  dd
})
tab_form <- bf(y ~ x + (1 | g)) + gaussian()

test_that("tibble input reproduces the data.frame fit (G2.7, G2.10)", {
  skip_if_not_installed("tibble")
  f_df <- frm(tab_form, data = tab_dd)
  f_tb <- frm(tab_form, data = tibble::as_tibble(tab_dd))
  expect_identical(as.numeric(logLik(f_tb)), as.numeric(logLik(f_df)))
  expect_identical(unlist(fixef(f_tb)), unlist(fixef(f_df)))
  p_df <- predict(f_df, newdata = tab_dd[1:5, ])
  p_tb <- predict(f_tb, newdata = tibble::as_tibble(tab_dd)[1:5, ])
  expect_identical(unname(p_tb), unname(p_df))
})

test_that("data.table input reproduces the data.frame fit (G2.7)", {
  skip_if_not_installed("data.table")
  f_df <- frm(tab_form, data = tab_dd)
  f_dt <- frm(tab_form, data = data.table::as.data.table(tab_dd))
  expect_identical(as.numeric(logLik(f_dt)), as.numeric(logLik(f_df)))
  expect_identical(unlist(fixef(f_dt)), unlist(fixef(f_df)))
})

test_that("columns with extra classes and attributes fit intact (G2.11)", {
  d <- tab_dd
  attr(d$x, "units") <- "mg"
  class(d$x) <- c("units_like", "numeric")
  f_cl <- frm(tab_form, data = d)
  f_df <- frm(tab_form, data = tab_dd)
  expect_identical(as.numeric(logLik(f_cl)), as.numeric(logLik(f_df)))
  expect_identical(unlist(fixef(f_cl)), unlist(fixef(f_df)))
})

test_that("list columns: unused ignored, used errors by name (G2.12)", {
  d <- tab_dd
  d$lc <- replicate(60, list(1:2), simplify = FALSE)
  f_lc <- frm(tab_form, data = d)
  f_df <- frm(tab_form, data = tab_dd)
  expect_identical(as.numeric(logLik(f_lc)), as.numeric(logLik(f_df)))
  expect_error(frm(bf(y ~ lc + (1 | g)) + gaussian(), data = d),
               "invalid type \\(list\\).*lc")
})

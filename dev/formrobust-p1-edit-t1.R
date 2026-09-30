# Punch round 1: tests for m5, m6, m7 (cens y2) and m8 in
# tests/testthat/test-aterm-expr.R. Record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/test-aterm-expr.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1('test_that("a term reads a constant from the formula environment", {
  k <- 2
  ae_same(frm(bf(y | weights(wt * k) ~ x), data = ae_data),
          frm(bf(y | weights(w2) ~ x), data = ae_data))
})',
'test_that("a term reads a constant from the formula environment", {
  k <- 2
  ae_same(frm(bf(y | weights(wt * k) ~ x), data = ae_data),
          frm(bf(y | weights(w2) ~ x), data = ae_data))
  # bare, as an expression and as a literal, the constant is one value:
  # brms refuses the first two (it reads the data alone), and frmtmb
  # reads all three alike, as the formula environment is read everywhere
  kt <- 12L
  d <- ae_data
  d$yt <- pmin(d$yb, 12L)
  f_lit <- frm(bf(yt | trials(12) ~ x), data = d, family = binomial())
  expect_identical(logLik(frm(bf(yt | trials(kt) ~ x), data = d,
                              family = binomial())), logLik(f_lit))
  expect_identical(logLik(frm(bf(yt | trials(kt + 0L) ~ x), data = d,
                              family = binomial())), logLik(f_lit))
})

test_that("a variable R finds only as a function is refused by name", {
  # without a column t, t is base::t() and the product failed with
  # "non-numeric argument to binary operator"
  expect_error(frm(bf(y | weights(t * 2) ~ x), data = ae_data),
               "R finds only the function t()", fixed = TRUE)
  # the guard-absent case: a column named t is read
  d <- ae_data
  d$t <- d$wt
  ae_same(frm(bf(y | weights(t * 2) ~ x), data = d),
          frm(bf(y | weights(w2) ~ x), data = d))
})

test_that("an expression that gives NA is refused, as brms refuses it", {
  expect_error(frm(bf(y | trunc(lb = ifelse(x > 1, NA, -5)) ~ x),
                   data = ae_data),
               "trunc(lb = ifelse(x > 1, NA, -5)) of response \'y\' is NA",
               fixed = TRUE)
  expect_error(frm(bf(yb | trials(ifelse(x > 1, NA, n + 1)) ~ x),
                   data = ae_data, family = binomial()),
               "is NA on", fixed = TRUE)
  expect_error(frm(bf(y | weights(ifelse(x > 1, NA, wt)) ~ x),
                   data = ae_data), "is NA on", fixed = TRUE)
  # the guard-absent case: an NA in the variable itself drops the row,
  # through na.action, as it does for any model variable
  d <- ae_data
  d$wt[1] <- NA
  f <- suppressMessages(frm(bf(y | weights(wt * 2) ~ x), data = d))
  expect_identical(nobs(f), nrow(d) - 1L)
})

test_that("the interval bound of cens() takes one value per row", {
  d <- ae_data
  d$cc <- ifelse(d$x > 1, "interval", "none")
  expect_error(frm(bf(y | cens(cc, max(y) + 10) ~ x), data = d),
               "the interval upper bound has 1 value(s) where the data have",
               fixed = TRUE)
  d$y2 <- max(d$y) + 10
  expect_identical(
    logLik(frm(bf(y | cens(cc, y2) ~ x), data = d)),
    logLik(frm(bf(y | cens(cc, y2 + 0) ~ x), data = d)))
})')
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\n")
close(con)

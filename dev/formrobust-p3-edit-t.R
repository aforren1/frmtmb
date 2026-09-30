# Round 3 (user decision, reversing m6): addition-term variables come
# from the data alone. Record of the edit of test-aterm-expr.R.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/test-aterm-expr.R"
x <- readLines(p)
i <- which(x == 'test_that("a term reads a constant from the formula environment", {')
j <- which(x == 'test_that("a variable R finds only as a function is refused by name", {')
stopifnot(length(i) == 1L, length(j) == 1L)
new <- c(
'test_that("a term reads its variables from the data alone, as brms", {',
'  # an outside value would be read again by predictions and refits, and',
'  # after saveRDS()/readRDS() in another session it is missing or',
'  # changed; on the response side, as trials(k), the change is silent',
'  k <- 2',
'  kt <- 12L',
'  d <- ae_data',
'  d$yt <- pmin(d$yb, 12L)',
'  msg <- "reads `k`, which is not a column of `data`"',
'  expect_error(frm(bf(y | weights(wt * k) ~ x), data = ae_data), msg,',
'               fixed = TRUE)',
'  expect_error(frm(bf(yc | rate(time * k) ~ x), data = ae_data,',
'                   family = poisson()), msg, fixed = TRUE)',
'  expect_error(frm(bf(yt | trials(kt) ~ x), data = d, family = binomial()),',
'               "reads `kt`, which is not a column of `data`", fixed = TRUE)',
'  expect_error(frm(bf(yt | trials(kt + 0L) ~ x), data = d,',
'                   family = binomial()),',
'               "Put `kt` in the data", fixed = TRUE)',
'  # the guard-absent cases: the same constant as a data column',
'  d$kt <- 12L',
'  f_lit <- frm(bf(yt | trials(12) ~ x), data = d, family = binomial())',
'  expect_identical(logLik(frm(bf(yt | trials(kt) ~ x), data = d,',
'                              family = binomial())), logLik(f_lit))',
'  expect_identical(logLik(frm(bf(yt | trials(kt + 0L) ~ x), data = d,',
'                              family = binomial())), logLik(f_lit))',
'  d$k <- 2',
'  ae_same(frm(bf(y | weights(wt * k) ~ x), data = d),',
'          frm(bf(y | weights(w2) ~ x), data = d))',
'  # function calls are not variables',
'  d$lw <- log(d$wt) + 1',
'  expect_identical(logLik(frm(bf(y | weights(exp(log(wt))) ~ x),',
'                              data = d)),',
'                   logLik(frm(bf(y | weights(wt) ~ x), data = d)))',
'  # and a constant in a predictor keeps the R convention',
'  kx <- 3',
'  d2 <- ae_data',
'  d2$x3 <- d2$x * 3',
'  expect_identical(logLik(frm(bf(y ~ I(x * kx)), data = ae_data)),',
'                   logLik(frm(bf(y ~ x3), data = d2)))',
'})',
'')
x <- c(x[seq_len(i - 1L)], new, x[j:length(x)])
writeLines(x, p)

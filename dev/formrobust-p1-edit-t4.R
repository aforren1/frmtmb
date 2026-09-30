# Punch round 1: update() tests (vacuous message check, call text,
# multivariate refusal). Record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/test-update-pool.R"
x <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1('  up2 <- NULL
  expect_message(
    up2 <- suppressMessages(update(f0, formula. = bf(y ~ a + b * x, a ~ 1,
                                                    nl = TRUE))),
    NA)',
'  msgs <- character()
  up2 <- withCallingHandlers(
    update(f0, formula. = bf(y ~ a + b * x, a ~ 1, nl = TRUE)),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(any(grepl("Replacing initial definitions of parameters a",
                        msgs, fixed = TRUE)))')
rep1('test_that("a plain formula keeps a linear model\'s dpar formulas", {
  f0 <- frm(bf(y ~ x, sigma ~ z), data = up_data)
  up <- update(f0, y ~ x + z)',
'test_that("a plain formula keeps a linear model\'s dpar formulas", {
  f0 <- frm(bf(y ~ x, sigma ~ z), data = up_data)
  up <- update(f0, y ~ x + z)
  # the call holds the pooled formula as the bf() call that builds it,
  # not the evaluated object with its family\'s closures
  expect_identical(deparse1(getCall(up)$formula),
                   paste0("frmtmb::bf(y ~ x + z, sigma ~ z, ",
                          "family = stats::gaussian(link = \\"identity\\"))"))
  expect_identical(logLik(eval(getCall(up))), logLik(up))')
x <- paste0(x, '

test_that("a complete formula on a multivariate fit is refused, as brms", {
  fm <- frm(bf(y ~ x, sigma ~ z) + bf(z ~ x), data = up_data)
  expect_error(update(fm, y ~ x + z),
               "Updating formulas of multivariate models is not yet",
               fixed = TRUE)
  expect_error(update(fm, mvbf(bf(y ~ x), bf(z ~ x))),
               "Updating formulas of multivariate models is not yet",
               fixed = TRUE)
  # the guard-absent case: the data can still be updated
  up <- update(fm, newdata = up_data[1:60, ])
  expect_identical(nobs(up), 60L)
})
')
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\n")
close(con)

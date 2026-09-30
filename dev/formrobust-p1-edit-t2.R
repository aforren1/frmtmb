# Punch round 1: bernoulli tests for m1 and m2. Record of the edit.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/"
edit <- function(f, old, new, eol = "\n") {
  p <- paste0(wt, f)
  x <- paste(readLines(p), collapse = "\n")
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <- sub(old, new, x, fixed = TRUE)
  con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = eol)
  close(con)
}
edit("test-bernoulli-coding.R",
'  # brms codes two fractions too; frmtmb refuses them as proportions
  # (a recorded divergence, lme4#682), and whole numbers pass above
  expect_error(cd(c(0.1, 0.9, 0.1)), "refuses fractions", fixed = TRUE)',
'  # two fractions are coded as brms codes them, with a warning, since
  # they are usually a proportion (lme4#682)
  pw <- "lie strictly between 0 and 1"
  allow_warnings(expect_identical(cd(c(0.1, 0.9, 0.1)), c(0, 1, 0)), pw,
                 require = pw)
  # the guard-absent cases: effect coding, and a value at 0 or 1, are
  # not a proportion and do not warn
  expect_no_warning(expect_identical(cd(c(-0.5, 0.5, 0.5)), c(0, 1, 1)))
  expect_no_warning(expect_identical(cd(c(0, 0.5, 0)), c(0, 1, 0)))')
edit("test-bernoulli-coding.R",
'test_that("brms\'s standata agrees", {
  skip_unless_brms()
  for (y in list(rep(-c(1, 2), 5), c(1, 2, 2, 1), c(5, 5, 5),',
'test_that("brms\'s standata agrees", {
  skip_unless_brms()
  for (y in list(rep(-c(1, 2), 5), c(1, 2, 2, 1), c(5, 5, 5),
                 c(-0.5, 0.5, 0.5), c(0.1, 0.9, 0.9),')
edit("test-bernoulli-coding.R",
'    fr <- frm(bf(y ~ 1), data = d, family = bernoulli(), dry_run = "frame")',
'    fr <- NULL
    allow_warnings(fr <- frm(bf(y ~ 1), data = d, family = bernoulli(),
                             dry_run = "frame"),
                   "lie strictly between 0 and 1")')
edit("test-bernoulli-coding.R",
'  # the coding travels with the fit\'s spec, which a refit reads',
'  # refit() takes the 0/1 codes simulate() returns; the response\'s own
  # values are refused by name, where they fitted a logLik of 3e304
  expect_error(refit(fm, bc_data$ym), "refitted on its 0/1 codes",
               fixed = TRUE)
  expect_identical(logLik(refit(fm, bc_data$y01)), logLik(ref))
  # the coding travels with the fit\'s spec, which a refit reads')
edit("test-open-issues.R",
'  expect_error(frm(bf(Y ~ 1 + (1 | g)) + bernoulli(), data = bd),
               "refuses fractions")',
'  pw <- "lie strictly between 0 and 1"
  allow_warnings(frm(bf(Y ~ 1 + (1 | g)) + bernoulli(), data = bd),
                 pw, require = pw)', eol = "\r\n")

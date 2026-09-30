# One-off edit of two tests that pinned refusals this lane lifts: the
# lower-case intercept (now brms's deprecated column of ones) and the
# bernoulli message (now brms's words). Kept as the record of the edit.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/"
edit <- function(f, old, new) {
  p <- paste0(wt, f)
  x <- paste(readLines(p), collapse = "\n")
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <- sub(old, new, x, fixed = TRUE)
  con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
  close(con)
}
edit("test-rsv-intercept.R",
  paste0('  expect_error(frm(bf(y ~ 0 + intercept + x), data = ri_data),\n',
         '               "deprecated spelling of the reserved variable `Intercept`")'),
  paste0('  # brms\'s deprecated lower-case spelling is read as brms reads it,\n',
         '  # a column of ones with a warning (test-brms-api-formrobust.R)\n',
         '  lw <- "Reserved variable name \'intercept\' is deprecated"\n',
         '  allow_warnings(\n',
         '    expect_identical(rownames(fixef(frm(bf(y ~ 0 + intercept + x),\n',
         '                                        data = ri_data))),\n',
         '                     c("intercept", "x")),\n',
         '    lw, require = lw)'))
edit("test-v14.R",
  '  expect_error(frm(bf(yc ~ x) + bernoulli(), data = dd), "0/1")',
  paste0('  # brms\'s words: two values of any kind are coded 0 and 1\n',
         '  expect_error(frm(bf(yc ~ x) + bernoulli(), data = dd),\n',
         '               "only two different values")'))

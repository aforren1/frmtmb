# One-off edit of two tests that pinned the bf()-time refusal of
# nl = TRUE without parameter formulas, which now happens at frm(), as
# in brms. Kept as the record of the edit.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/"
edit <- function(f, old, new) {
  p <- paste0(wt, f)
  x <- paste(readLines(p), collapse = "\n")
  stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
  x <- sub(old, new, x, fixed = TRUE)
  con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
  close(con)
}
edit("test-nl.R",
  '  expect_error(bf(y ~ a * exp(-b * x), nl = TRUE), "parameter formula")',
  paste0('  # refused when the model is assembled, as brms refuses it at brm():\n',
         '  # the formulas may still arrive with + lf() (lane formrobust)\n',
         '  expect_error(frm(bf(y ~ a * exp(-b * x), nl = TRUE), data = NULL,\n',
         '                   dry_run = "spec"),\n',
         '               "nonlinear-parameter formula")'))
edit("test-parse.R",
  '  expect_error(bf(y ~ a * exp(b * x), nl = TRUE), "parameter formula")',
  paste0('  # nl = TRUE without a parameter formula waits for frm(), as brms\n',
         '  # waits for brm()\n',
         '  expect_error(frm(bf(y ~ a * exp(b * x), nl = TRUE), data = NULL,\n',
         '                   dry_run = "spec"),\n',
         '               "nonlinear-parameter formula")'))

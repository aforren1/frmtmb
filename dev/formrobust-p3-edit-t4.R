# Round 3: the absent case of the se() refusal now meets the addition-
# term data rule first. Record of the edit.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/test-brms-parity-defects.R"
x <- paste(readLines(p), collapse = "\n")
old <- '  # the absent case: a family that reads se() gets the missing column
  expect_error(frm(y | se(sei) ~ x, data = dd, family = gaussian()),
               "The model uses `sei`", fixed = TRUE)'
new <- '  # the absent case: a family that reads se() gets the missing column,
  # named by the addition-term rule (a variable must be in the data)
  expect_error(frm(y | se(sei) ~ x, data = dd, family = gaussian()),
               "reads `sei`, which is not a column of `data`", fixed = TRUE)'
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)

# One-off for the fix round: replace the runaway test of
# test-quadrature-defects.R (from its test_that line up to the next
# test_that) with dev/quadseam-fr-newtests.txt.
f <- "C:/Users/adf44/source/r/frmtmb-wt-quadseam/tests/testthat/test-quadrature-defects.R"
x <- readLines(f)
a <- which(x == paste0('test_that("quadrature does not stop a runaway at ',
                       'a false maximum", {'))
b <- which(x == paste0('test_that("a quadrature sd at zero is a boundary ',
                       'fit, as under Laplace", {'))
stopifnot(length(a) == 1L, length(b) == 1L, b > a)
new <- readLines(
  "C:/Users/adf44/source/r/frmtmb-wt-quadseam/dev/quadseam-fr-newtests.txt")
x <- c(x[seq_len(a - 1L)], new, "", x[b:length(x)])
writeLines(x, f)
cat("spliced", length(new), "lines\n")

# expect_vector_equal(), copied from frmtmb's tests/testthat/helper-
# reference.R for the tests that moved here with car(): a helper the
# suites share is cheaper to copy than to promise as API.
expect_vector_equal <- function(x, y, tol) {
  expect_equal(length(x), length(y))
  expect_lt(max(abs(unname(x) - unname(y))), tol)
}

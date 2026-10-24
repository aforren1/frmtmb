# Copied from frmtmb's helper-reference.R rather than reached for, as
# frmtmb.sample's helper-sampling.R does: a three-line expectation is
# cheaper to duplicate than to make into API. The cox() tests that
# moved here from frmtmb use it.
expect_vector_equal <- function(x, y, tol) {
  expect_equal(length(x), length(y))
  expect_lt(max(abs(unname(x) - unname(y))), tol)
}

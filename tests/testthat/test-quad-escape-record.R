# A quadrature fit escapes a family's stationary point on the Laplace
# objective (quad_fit()) and again on the rule (fit_assembled()); the
# second escape's record replaced the first's, so a fit the Laplace escape
# moved and the rule's left alone reported a gain of 0
# (stationary_escape_keep(), release review of 0.82.0, minor 6).

test_that("the escape that moved a quadrature fit keeps its record", {
  lap <- c(starts = 2, gain = 24.27)
  # the rule's escape ran and gained nothing: the Laplace record stays
  opt <- list(objective = 1, stationary_escape = c(starts = 2, gain = 0))
  got <- stationary_escape_keep(opt, lap)
  expect_identical(got$stationary_escape, lap)
  expect_identical(got$stationary_escape_laplace, lap)
  # the rule's escape moved the fit: its record stays, the Laplace one
  # beside it
  rule <- c(starts = 2, gain = 0.5)
  opt$stationary_escape <- rule
  got <- stationary_escape_keep(opt, lap)
  expect_identical(got$stationary_escape, rule)
  expect_identical(got$stationary_escape_laplace, lap)
})

# gp_chol_pd(), the positive-definiteness check of the inner Hessian that
# gates the second reading of a gp() body (R/gp-start.R). CHOLMOD factors
# a matrix with a NaN or an infinite entry without a word, so the check
# reads the entries first; on a finite matrix only CHOLMOD's "not
# positive" warning or error is a FALSE verdict.

ipd_matrix <- function(d, off = 0.5) {
  m <- diag(d)
  m[1, 2] <- m[2, 1] <- off
  Matrix::forceSymmetric(methods::as(Matrix::Matrix(m, sparse = TRUE),
                                     "CsparseMatrix"))
}

test_that("a non-finite inner Hessian is not positive definite", {
  expect_false(gp_chol_pd(ipd_matrix(c(NaN, 2, 3))))
  expect_false(gp_chol_pd(ipd_matrix(c(2, Inf, 3))))
  expect_false(gp_chol_pd(ipd_matrix(c(2, 2, 3), off = NaN)))
})

test_that("a finite matrix is judged by its factor, in silence", {
  expect_true(gp_chol_pd(ipd_matrix(c(2, 2, 3))))
  # the leading 2 x 2 block has determinant 1 - 4 < 0
  expect_no_warning(got <- gp_chol_pd(ipd_matrix(c(1, 1, 3), off = 2)))
  expect_false(got)
  expect_no_warning(got <- gp_chol_pd(ipd_matrix(c(-1, 2, 3))))
  expect_false(got)
})

test_that("another warning inside the factor reaches the caller", {
  local_mocked_bindings(
    Cholesky = function(...) {
      warning("a deprecation", call. = FALSE)
      invisible(NULL)
    },
    .package = "Matrix")
  got <- allow_warnings(gp_chol_pd(ipd_matrix(c(2, 2, 3))), "a deprecation",
                        require = "a deprecation")
  expect_true(got)
})

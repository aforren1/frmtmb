# Fix round 1 edits to frmtmb.latent's test-mixture-mvn.R; each old
# string must occur exactly once.
f <- "extensions/frmtmb.latent/tests/testthat/test-mixture-mvn.R"
x <- gsub("\r\n", "\n", readChar(f, file.info(f)$size, useBytes = TRUE),
          fixed = TRUE)
rep1 <- function(old, new) {
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(n, " matches for: ", old)
  x <<- sub(old, new, x, fixed = TRUE)
}
rep1("# in its own tests. Three frmtmb internals the old file reached are
# replaced by the local helpers below: make_start(), us_theta_cor() and
# the expect_vector_equal() of frmtmb's helper-reference.R.",
"# in its own tests. Three frmtmb internals the old file reached are
# replaced: make_start() and the expect_vector_equal() of frmtmb's
# helper-reference.R by the local helpers below, and us_theta_cor() by
# the exported frmtmb::us_theta().")
rep1("# The correlation segment of a us() parameter vector from a correlation
# matrix: the inverse of the map frmtmb::us_sigma() documents. C = Lc
# Lc' with Lc lower-triangular has unit-norm rows, so dividing row i by
# Lc[i, i] recovers the unit-diagonal factor.
mvn_theta_cor <- function(C) {
  Lc <- t(chol(C))
  L <- Lc / diag(Lc)
  L[lower.tri(L)]
}

",
"")
rep1("  us <- function(S) {
    c(log(sqrt(diag(S))), mvn_theta_cor(stats::cov2cor(S)))
  }",
"  us <- function(S) us_theta(S)")
rep1("test_that(\"method = 'predict' refuses mixture_mvn(), which draws the response whole\", {",
"test_that(\"method = 'predict' refuses a draw of the whole response\", {")
writeChar(x, f, eos = NULL, useBytes = TRUE)
cat("done\n")

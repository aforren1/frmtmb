p <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.spatial/tests/testthat/test-sar.R"
txt <- paste(readLines(p), collapse = "\n")
rep1 <- function(old, new) {
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop("matches: ", n, " for ", substr(old, 1, 60))
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1("  # and a maximum
  X <- cbind(1, d$x)
  ev <- eigen(W, only.values = TRUE)$values
  lo <- 1 / min(Re(ev))
  hi <- 1 / max(Re(ev))
  nll <- function(p) {
    r <- lo + (hi - lo) / (1 + exp(-p[4]))
    Ap <- diag(s$n) - r * W
    -dense_mvn(d$ylag, solve(Ap, X %*% p[1:2]),
               exp(2 * p[3]) * solve(crossprod(Ap)))
  }",
"  # and a maximum: the same density in its change-of-variables form,
  # log|A| from a dense determinant, optimized on its own
  X <- cbind(1, d$x)
  ev <- eigen(W, only.values = TRUE)$values
  lo <- 1 / min(Re(ev))
  hi <- 1 / max(Re(ev))
  nll <- function(p) {
    r <- lo + (hi - lo) / (1 + exp(-p[4]))
    Ap <- diag(s$n) - r * W
    e <- as.vector(Ap %*% d$ylag - X %*% p[1:2])
    -(sum(stats::dnorm(e, 0, exp(p[3]), log = TRUE)) +
        as.numeric(determinant(Ap)$modulus))
  }")
rep1("  expect_equal(-fit$obj$fn(p), ref, tolerance = 1e-10)",
     "  expect_equal(-fit$obj$fn(p), unname(ref), tolerance = 1e-10)")
rep1("  d$ylag <- d$ylag + rnorm(7, 0, 0.5)[d$g]",
     "  d$ylag <- d$ylag + rnorm(7, 0, 1.5)[d$g]")
writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)

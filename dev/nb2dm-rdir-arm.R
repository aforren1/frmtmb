# test-simplex-rdirichlet.R on the lane build with the probability-scale
# Dirichlet draw of rellib-r10 put back into simplex_rdirichlet(), and
# the row-failure counts of that draw at phi 0.01 and 1e-3.
#   Rscript dev/nb2dm-rdir-arm.R > dev/nb2dm-out/rdir-arm.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb)})
old <- function(alpha, cats) {
  G <- matrix(stats::rgamma(length(alpha), shape = alpha), nrow(alpha))
  out <- G / rowSums(G)
  colnames(out) <- cats
  out
}
new <- frmtmb:::simplex_rdirichlet
p <- exp(c(0, 0.3, -0.2)) / sum(exp(c(0, 0.3, -0.2)))
n <- 20000
for (phi in c(0.01, 1e-3)) {
  a <- matrix(rep(p * phi, each = n), n)
  set.seed(5)
  b_old <- sum(!is.finite(rowSums(old(a, letters[1:3]))))
  set.seed(5)
  b_new <- sum(!is.finite(rowSums(new(a, letters[1:3]))))
  cat(sprintf("phi %g: NaN rows of %d, probability-scale draw %d, log-scale draw %d\n",
              phi, n, b_old, b_new))
}
utils::assignInNamespace("simplex_rdirichlet", old, "frmtmb")
r <- as.data.frame(test_file("tests/testthat/test-simplex-rdirichlet.R",
                             package = "frmtmb", env = test_env("frmtmb"),
                             reporter = "silent", stop_on_failure = FALSE))
print(r[, c("test", "failed", "error")])

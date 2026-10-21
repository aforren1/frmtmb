# Where are the Gamma CDF pair's third derivatives not finite? Seed 7
# draws as dev/cdfs-grad.R makes them, then RTMB::pgamma() alone.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib", "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
ns <- asNamespace("frmtmb")
"c" <- RTMB::ADoverload("c")
q <- c(0.01, 0.3, 1, 2.5, 6, 20)
set.seed(7)
pts <- cbind(runif(200, -2, 3), runif(200, -2, 6))
for (j in seq_along(q)) {
  f <- function(p) { lp <- ns$gamma_cdf_log(q[j], exp(p[1]), exp(p[2])); c(lp$lF, lp$lS) }
  T3 <- MakeTape(f, c(log(2), log(3)))$jacfun()$jacfun()
  bad <- which(apply(pts, 1, function(p) !all(is.finite(T3$jacobian(p)))))
  if (length(bad)) {
    x <- q[j] * exp(pts[bad, 2]) / exp(pts[bad, 1]); a <- exp(pts[bad, 2])
    cat("q", q[j], "bad", length(bad), "\n"); print(cbind(a = a, x = x, x_over_ap1 = x / (a + 1), d = (x - a - 1) / sqrt(a)))
  }
}
g <- function(p) log(RTMB::pgamma(exp(p[1]), shape = exp(p[2])))
T3 <- MakeTape(g, c(0, 0))$jacfun()$jacfun()
xs <- cbind(runif(2000, -3, 4), runif(2000, -2, 6))
bad <- apply(xs, 1, function(p) !all(is.finite(T3$jacobian(p))))
cat("plain log(RTMB::pgamma) third derivative non-finite:", sum(bad), "of", nrow(xs), "\n")
print(head(cbind(x = exp(xs[bad, 1]), a = exp(xs[bad, 2])), 10))

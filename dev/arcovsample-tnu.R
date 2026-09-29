# Lane wt-arcovsample: which side of the student rescor disagreement is
# wrong, and at what nu. The objective uses mvt_std_loglik(), which goes
# through lgamma_shift_diff(); rescor_row_loglik() writes
# lgamma((nu + K) / 2) - lgamma(nu / 2) and log(nu * pi) as they read.
#   Rscript dev/arcovsample-tnu.R
.libPaths(c("C:/Users/adf44/source/r/wt-arcovsample-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))

set.seed(7)
K <- 2L
n <- 5L
Z <- matrix(stats::rnorm(n * K), n, K)
C <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
q <- rowSums((Z %*% solve(C)) * Z)
ldet <- as.numeric(determinant(C, logarithm = TRUE)$modulus)

naive <- function(nu) {
  lgamma((nu + K) / 2) - lgamma(nu / 2) - K / 2 * log(nu * pi) -
    ldet / 2 - (nu + K) / 2 * log1p(q / nu)
}
stable <- function(nu) {
  frmtmb:::lgamma_shift_diff(nu / 2, K / 2) -
    K / 2 * (log(nu) + log(pi)) - ldet / 2 -
    (nu + K) / 2 * log1p(q / nu)
}
# the gaussian limit, which both must approach as nu grows
gaussian <- function() {
  -K / 2 * log(2 * pi) - ldet / 2 - q / 2
}
# the objective's own function, summed over the rows
obj <- function(nu) frmtmb:::mvt_std_loglik(Z, C, nu)

old <- options(digits = 12)
cat(sprintf("%12s %18s %18s %18s %14s\n", "nu", "sum naive",
            "sum stable", "objective", "naive-obj"))
for (nu in c(5, 30, 1e3, 1e5, 1e7, 1e10, 1e20, 1e100, 1e250, 1e300,
             1e305, 1e306)) {
  a <- sum(naive(nu))
  b <- sum(stable(nu))
  o <- as.numeric(obj(nu))
  cat(sprintf("%12.0e %18.8f %18.8f %18.8f %14.3e\n", nu, a, b, o, a - o))
}
cat("\ngaussian limit (sum): ", format(sum(gaussian())), "\n")
cat("stable at nu = 1e306: ", format(sum(stable(1e306))), "\n")
cat("naive  at nu = 1e306: ", format(sum(naive(1e306))), "\n")
options(old)

# Reviewer: the small-phi failure of dev/nb2dm-rev-smallphi.R, counted.
# How often does the family's simulator return a NaN row as phi falls,
# and does a log-scale Dirichlet draw (log G = log Gamma(a + 1) +
# log(U) / a) remove it? 20000 rows per phi, K = 3, trials 12, seed 5.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
P <- exp(c(0, 0.3, -0.2)) / sum(exp(c(0, 0.3, -0.2)))
naive_bad <- function(alpha) {
  G <- matrix(stats::rgamma(length(alpha), shape = alpha), nrow(alpha))
  sum(!is.finite(rowSums(G / rowSums(G))))
}
log_draw <- function(alpha) {
  lg <- log(matrix(stats::rgamma(length(alpha), shape = alpha + 1),
                   nrow(alpha))) +
    log(matrix(stats::runif(length(alpha)), nrow(alpha))) / alpha
  e <- exp(lg - apply(lg, 1, max))
  e / rowSums(e)
}
n <- 20000
for (phi in c(1, 0.1, 0.05, 0.02, 0.01, 1e-3)) {
  set.seed(5)
  alpha <- matrix(rep(P * phi, each = n), n)
  bad <- naive_bad(alpha)
  set.seed(5)
  pr <- log_draw(alpha)
  # at small phi each row's probabilities sit at a vertex with
  # probability p_k: the share of rows whose largest p exceeds 0.99
  cat(sprintf(paste("phi %-6g NaN rows, rgamma draw: %5d of %d;",
                    "log-scale draw: %d; vertex share %.3f; mean p",
                    "%.3f %.3f %.3f (target %.3f %.3f %.3f)\n"),
              phi, bad, n, sum(!is.finite(rowSums(pr))),
              mean(apply(pr, 1, max) > 0.99), colMeans(pr)[1],
              colMeans(pr)[2], colMeans(pr)[3], P[1], P[2], P[3]))
}
fam <- dirichlet_multinomial(K = 3)
probe("the family's simulator at phi 0.02, 20000 rows", {
  set.seed(6)
  x <- fam$sim(list(mu2 = rep(0.3, n), mu3 = rep(-0.2, n),
                    phi = rep(0.02, n)), list(trials = rep(12, n)), n)
  dim(x)
})

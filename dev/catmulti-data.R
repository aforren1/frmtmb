# Data generators shared by lane catmulti's scripts (brms fits and the
# frmtmb comparisons read the same rows).

# A three-level factor `f` (levels a, b, c) from a multinomial logit on x.
catmulti_cat_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$f <- factor(apply(P, 1, function(p) {
    sample(c("a", "b", "c"), 1, prob = p)
  }), levels = c("a", "b", "c"))
  dd
}

# A three-column count matrix `Yc` (y1, y2, y3) with N = 10 trials a row.
catmulti_mult_data <- function(seed, n = 200) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -0.4 - 0.5 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  Yc <- t(vapply(seq_len(n), function(i) {
    as.numeric(stats::rmultinom(1, 10, P[i, ]))
  }, numeric(3)))
  colnames(Yc) <- c("y1", "y2", "y3")
  dd$Yc <- Yc
  dd$N <- 10L
  dd
}

# A composition `Y` (a, b, c) and a gaussian `w`, for mvbf().
catmulti_mv_data <- function(seed, n = 300) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.4 + 0.5 * dd$x, -0.3 - 0.4 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  G <- matrix(stats::rgamma(n * 3, shape = 6 * P), n, 3)
  Y <- G / rowSums(G)
  colnames(Y) <- c("a", "b", "c")
  dd$Y <- Y
  dd$w <- 1 + 0.5 * dd$x + stats::rnorm(n)
  dd
}

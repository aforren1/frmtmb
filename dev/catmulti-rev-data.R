# Reviewer data for lane catmulti: categorical responses with a grouping
# factor and per-category correlated intercepts, and a multinomial with
# varying trials.

# K categories, letters[1:K]; category 1 the generating reference.
rev_cat_data <- function(seed, K, n = 400, ng = 25) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(sample(paste0("g", seq_len(ng)), n, TRUE)))
  L <- t(chol(0.5 * diag(K - 1L) + 0.3))
  u <- matrix(stats::rnorm(ng * (K - 1L)), ng) %*% t(L)
  gi <- as.integer(dd$g)
  b0 <- seq(0.3, -0.4, length.out = K - 1L)
  b1 <- seq(0.8, -0.6, length.out = K - 1L)
  E <- cbind(0, sapply(seq_len(K - 1L), function(k) {
    b0[k] + b1[k] * dd$x + u[gi, k]
  }))
  P <- exp(E) / rowSums(exp(E))
  dd$y <- factor(apply(P, 1, function(p) sample(letters[1:K], 1, prob = p)),
                 levels = letters[1:K])
  dd
}

# Three categories y1..y3, trials N varying 1..30 a row.
rev_mult_data <- function(seed, n = 250) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  E <- cbind(0, 0.3 + 0.8 * dd$x, -1.5 - 1.2 * dd$x)
  P <- exp(E) / rowSums(exp(E))
  dd$N <- sample(1:30, n, TRUE)
  Yc <- t(vapply(seq_len(n), function(i) {
    as.numeric(stats::rmultinom(1, dd$N[i], P[i, ]))
  }, numeric(3)))
  colnames(Yc) <- c("y1", "y2", "y3")
  dd$Yc <- Yc
  dd$y1 <- Yc[, 1]; dd$y2 <- Yc[, 2]; dd$y3 <- Yc[, 3]
  dd
}

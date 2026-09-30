# Data shared by test-mo-terms.R and the brms log-density tier's
# ordinal-with-mo() row: one monotonic variable in two terms whose
# shapes a single simplex cannot hold at once.

mo_terms_data <- function(seed = 21, n = 400) {
  set.seed(seed)
  inc <- sample(0:3, n, replace = TRUE)
  w <- sample(0:2, n, replace = TRUE)
  z <- rnorm(n)
  # the main effect rises early, the interaction late: two shapes a
  # single simplex cannot hold at once
  main <- c(0, 0.7, 0.9, 1)[inc + 1]
  slope <- c(0, 0.05, 0.15, 1)[inc + 1]
  data.frame(y = 1 + 2 * main + (0.3 + 1.5 * slope) * z + rnorm(n, 0, 0.6),
             inc = inc, w = w, z = z, z2 = rnorm(n))
}


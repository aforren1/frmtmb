f <- "R/compat.R"
x <- readLines(f)
old1 <- 'Use type = \\"pearson\\", which divides by the marginal residual SD (the diagonal of the residual covariance is sigma^2, because R is unit-diagonal), or dharma_residuals().'
new1 <- 'Use type = \\"pearson\\", which divides by the marginal residual SD, or dharma_residuals().'
old2 <- '\\"response\\" and \\"pearson\\" are unchanged, and pearson is still the right standardization: R is unit-diagonal, so the marginal residual SD is sigma. They do NOT'
new2 <- '\\"response\\" is y - mu, and \\"pearson\\" divides it by the marginal residual SD of the row: sigma under cosy() and unstr(), and sigma times the ratio of the process sd to the innovation sd under ar(), ma() and arma(), whose sigma is the innovation sd, as in brms. They do NOT'
stopifnot(sum(grepl(old1, x, fixed = TRUE)) == 1, sum(grepl(old2, x, fixed = TRUE)) == 1)
x <- sub(old1, new1, x, fixed = TRUE)
x <- sub(old2, new2, x, fixed = TRUE)
writeLines(x, f)
f <- "R/autocor-loglik.R"
x <- readLines(f)
i <- grep("^#' form `e'Qe` and group size `k`. With `D = diag", x)
x <- c(x[seq_len(i - 1)],
  "#' form `e'Qe` and group size `k`. With `D = diag(sigma)`, `sigma` the",
  "#' marginal sd of each row, and `R` the pattern's correlation,",
  "#' `Q = D^-1 R^-1 D^-1`, so on the standardized residuals `z` the",
  "#' pattern needs one `k x k` inverse for all its groups.",
  x[(i + 4):length(x)])
writeLines(x, f)

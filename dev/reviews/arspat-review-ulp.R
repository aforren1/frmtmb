# Review of lane arspat: the mi() test's cross-fit identity (logLik of the
# mi() fit against the fit on the dropped rows, 1e-10 relative) under a
# few-ulp perturbation of y, 20 draws, gaussian and student.
source("dev/reviews/arspat-review-common.R")
acse_data <- function(seed = 3, G = 10, K = 5, ar = 0.5, tdf = NULL) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  n <- nrow(d)
  d$x <- stats::rnorm(n)
  d$s <- stats::runif(n, 0.2, 0.6)
  e <- as.vector(vapply(seq_len(G), function(i) {
    as.vector(stats::arima.sim(list(ar = ar), K, n.start = 100, sd = 0.7))
  }, numeric(K)))
  if (!is.null(tdf)) e <- e * rep(sqrt(tdf / stats::rchisq(G, tdf)), each = K)
  d$y <- 1 + 0.5 * d$x + e + stats::rnorm(n, 0, d$s)
  d
}
for (fam in list(gaussian(), student())) {
  d0 <- acse_data(seed = 6, tdf = if (fam$family == "student") 5)
  miss <- c(3, 12, 26, 49)
  worst <- 0
  for (k in 0:19) {
    d <- d0
    d$y <- d$y * (1 + 2e-15 * k)
    d$y[miss] <- NA
    fit <- suppressMessages(frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)),
                                family = fam, data = d))
    f0 <- suppressMessages(frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = fam,
                               data = d[-miss, ]))
    r <- abs(as.numeric(logLik(fit)) - as.numeric(logLik(f0))) /
      abs(as.numeric(logLik(f0)))
    worst <- max(worst, r)
  }
  cat(sprintf("%s: worst relative logLik difference over 20 draws %.3e (bound 1e-10)\n",
              fam$family, worst))
}

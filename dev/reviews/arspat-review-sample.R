# Review of lane arspat: frm_sample() on ar(cov = TRUE) + mi() and on
# se(): every missing value is sampled, and the draws' posterior mean of
# the coefficients sits near the ML fit.
source("dev/reviews/arspat-review-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
set.seed(23)
G <- 10; K <- 6
d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
n <- nrow(d)
d$x <- rnorm(n); d$s <- runif(n, 0.2, 0.5)
e <- as.vector(vapply(seq_len(G), function(i)
  as.vector(arima.sim(list(ar = 0.6), K, sd = 0.8)), numeric(K)))
d$y <- 1 + 0.5 * d$x + e
miss <- c(3, 4, 14, 31, 36, 55)
d$y[miss] <- NA
chk("frm_sample ar + mi", {
  f <- frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), data = d)
  s <- frm_sample(f, chains = 2, iter = 1000, refresh = 0, seed = 3)
  v <- variables(s)
  cat("variables:", paste(v, collapse = " "), "\n")
  cat("Ymi-like variables:", sum(grepl("mi|Ymi", v)), "of", length(miss), "\n")
  m <- as.matrix(s)
  cat("posterior mean b_x", mean(m[, grep("^b_x$", colnames(m))]),
      " ML", fixef(f)[2, 1], "\n")
  NULL
})
chk("frm_sample ar + se, log_lik and loo", {
  d2 <- d; d2$y[miss] <- 1 + 0.5 * d2$x[miss]
  f <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)), data = d2)
  s <- frm_sample(f, chains = 2, iter = 1000, refresh = 0, seed = 3)
  ll <- log_lik(s)
  cat("log_lik dim", dim(ll), " finite", all(is.finite(ll)), "\n")
  pp <- posterior_predict(s, ndraws = 200)
  cat("posterior_predict dim", dim(pp), "\n")
  NULL
})

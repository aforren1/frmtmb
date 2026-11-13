# Review of lane arspat: predict() Est.Error under se() with the parameter
# draws off (propagate_error = FALSE), and 4320 against a hand
# computation of brms's cov = FALSE one-step mean.
source("dev/reviews/arspat-review-common.R")
set.seed(7)
G <- 40; K <- 5
d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
n <- nrow(d)
d$x <- rnorm(n); d$x2 <- rnorm(n); d$s <- runif(n, 0.2, 0.9)
e <- as.vector(vapply(seq_len(G), function(i)
  as.vector(arima.sim(list(ar = 0.6), K, sd = 0.6)), numeric(K)))
d$y <- 1 + 0.6 * d$x + e * exp(0.3 * d$x2) + rnorm(n, 0, d$s)
chk("predict Est.Error under se(), propagate_error = FALSE, 20000 draws", {
  f <- frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE),
              sigma ~ x2), data = d)
  p <- predict(f, ndraws = 20000, propagate_error = FALSE)
  dp <- ns$eval_dpars(f)[[1L]]
  phi <- unname(ac_nat(f)[1])
  v <- dp$sigma^2 / (1 - phi^2) + d$s^2
  r <- p[, "Estimate"] * 0 + p[, "Est.Error"] / sqrt(v)
  cat("ratio Est.Error / sqrt(S_ii): range", format(range(r), digits = 4),
      " mean", format(mean(r), digits = 5), " (MC sd per row ~",
      format(sqrt(1 / 40000), digits = 2), ")\n")
  NULL
})
chk("4320 against the hand one-step mean", {
  set.seed(2)
  dd <- data.frame(t = rep(1:20, 3), h = factor(rep(1:3, each = 20)),
                   x = rnorm(60))
  dd$y <- dd$x + rnorm(60)
  f0 <- frm(bf(y ~ x + ar(time = t, gr = h)), data = dd)
  nd <- dd[c(1, 2, 3, 21, 22), ]
  p <- marginaleffects::predictions(f0, newdata = nd)
  b <- fixef(f0)[, 1]
  phi <- unname(ac_nat(f0)[1])
  mu <- b[1] + b[2] * nd$x
  hand <- mu
  hand[2] <- mu[2] + phi * (nd$y[1] - mu[1])
  hand[3] <- mu[3] + phi * (nd$y[2] - mu[2])
  hand[5] <- mu[5] + phi * (nd$y[4] - mu[4])
  cat("predictions:", format(p$estimate, digits = 10), "\n")
  cat("hand       :", format(unname(hand), digits = 10), "\n")
  cat(sprintf("max abs diff %.3e\n", max(abs(p$estimate - hand))))
  # newdata with no padding: unchanged route
  nd2 <- dd[dd$h == "1", ]
  p2 <- marginaleffects::predictions(f0, newdata = nd2)
  f2 <- fitted(f0, newdata = nd2)[, "Estimate"]
  cat(sprintf("no-padding newdata max abs diff to fitted: %.3e (n %d)\n",
              max(abs(p2$estimate - f2)), nrow(p2)))
  # slopes and comparisons on a cov = FALSE fit
  print(marginaleffects::avg_slopes(f0, variables = "x"))
  NULL
})

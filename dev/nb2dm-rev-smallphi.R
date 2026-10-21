# Reviewer: a dirichlet_multinomial fit whose phi is small (seed 41,
# data drawn at phi = 0.02, n = 200, K = 3): does simulate() / predict()
# work, since the Dirichlet draw underflows at small concentrations?
# Then hmm() with a negbinomial2 emission against negbinomial.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
set.seed(41)
n <- 200
d <- data.frame(x = rnorm(n))
P <- exp(cbind(0, 0.3 + 0.4 * d$x, -0.2))
P <- P / rowSums(P)
d$N <- sample(5:30, n, TRUE)
d$Y <- t(vapply(seq_len(n), function(i) {
  # a Dirichlet draw on the log scale, so it does not underflow
  lg <- log(stats::rgamma(3, P[i, ] * 0.02 + 1)) +
    log(stats::runif(3)) / (P[i, ] * 0.02)
  p <- exp(lg - max(lg))
  as.numeric(rmultinom(1, d$N[i], p / sum(p)))
}, numeric(3)))
colnames(d$Y) <- c("a", "b", "c")
fit <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = d)
probe("phi estimate", exp(fixef_by_dpar(fit)$phi))
probe("simulate()", {
  s <- simulate(fit, nsim = 5, seed = 1)
  vapply(s, function(m) all(rowSums(m) == d$N), NA)
})
probe("predict()", {
  set.seed(2)
  dim(predict(fit, ndraws = 500))
})
suppressMessages(library(frmtmb.latent))
probe("hmm() negbinomial2 vs negbinomial emissions", {
  set.seed(12)
  dd <- data.frame(id = 1, t = 1:300)
  s <- integer(300)
  s[1] <- 1L
  G <- matrix(c(0.92, 0.08, 0.15, 0.85), 2, 2, byrow = TRUE)
  for (t in 2:300) s[t] <- sample.int(2, 1, prob = G[s[t - 1], ])
  dd$y <- rnbinom(300, mu = c(2, 9)[s], size = 3)
  h1 <- frm(bf(y ~ 1), family = hmm(K = 2, negbinomial(), time = t,
                                    group = id), data = dd)
  h2 <- frm(bf(y ~ 1), family = hmm(K = 2, brmsfamily("negbinomial2"),
                                    time = t, group = id), data = dd)
  c(nb = logLik(h1), nb2 = logLik(h2))
})

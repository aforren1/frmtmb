# Reviewer: odds and ends. dirichlet_multinomial's simulator at a tiny
# phi (brms's rdirichlet underflows there); frmtmb.latent hmm() with a
# negbinomial2 emission against negbinomial; family = "<name>" strings.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
fam <- dirichlet_multinomial(K = 3)
for (phi in c(1e-6, 1e-3, 0.05)) {
  probe(paste("sim at phi", phi), {
    set.seed(1)
    n <- 2000
    x <- fam$sim(list(mu2 = rep(0.3, n), mu3 = rep(-0.2, n),
                      phi = rep(phi, n)), list(trials = rep(8, n)), n)
    p <- exp(c(0, 0.3, -0.2)) / sum(exp(c(0, 0.3, -0.2)))
    # at a tiny phi a row puts all its trials in one category, with
    # probability p_k
    c(any_na = anyNA(x), rows_ok = all(rowSums(x) == 8),
      share_all_in_one = mean(apply(x, 1, max) == 8),
      p_all_in_a = mean(x[, 1] == 8), p_a = p[1])
  })
}
probe("brms rdirichlet at alpha 1e-6 * p", {
  set.seed(1)
  a <- matrix(rep(1e-6 * c(0.3, 0.4, 0.3), each = 2000), 2000)
  r <- brms::rdirichlet(2000, a[1, ])
  c(nan_rows = sum(!is.finite(rowSums(r))))
})
probe("family = 'negbinomial2' and 'dirichlet_multinomial' strings", {
  set.seed(2)
  d <- data.frame(x = rnorm(100))
  d$y <- rnbinom(100, mu = 3, size = 2)
  d$N <- 10L
  d$Y <- t(sapply(1:100, function(i) rmultinom(1, 10, c(.2, .3, .5))))
  c(logLik(frm(y ~ x, family = "negbinomial2", data = d)),
    logLik(frm(Y | trials(N) ~ x, family = "dirichlet_multinomial", data = d)))
})
suppressMessages(library(frmtmb.latent))
probe("hmm() negbinomial2 vs negbinomial emissions", {
  set.seed(3)
  n <- 400
  st <- rep(rep(1:2, each = 20), length.out = n)
  d <- data.frame(t = seq_len(n), id = 1,
                  y = rnbinom(n, mu = c(2, 9)[st], size = 3))
  h1 <- hmm(y ~ 1, data = d, nstates = 2, family = negbinomial(),
            time = "t", id = "id")
  h2 <- hmm(y ~ 1, data = d, nstates = 2,
            family = brmsfamily("negbinomial2"), time = "t", id = "id")
  c(nb = logLik(h1), nb2 = logLik(h2))
})

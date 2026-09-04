# Helpers the sampling suite shares.
#
# `sampler_gates_on()` came from frmtmb's helper-reference.R with the
# whole of its user base: every one of the twelve files that consulted
# it was a sampling file, so it left rather than being duplicated. The
# comment below is its original one, unchanged, because the reasons
# have not changed.

# Chain-agreement gates: assertions that compare NUTS chain output to a
# reference quantity (the ML fit, a Wald SE, another chain). Two
# reasons they carry a switch, learned in that order:
# (1) a seeded Stan chain is not platform-deterministic (a 1.11-sd
# mean gap was measured on macOS CI), so agreement bands must not
# assume this machine's chains; and
# (2) what looked like worse chain luck in pkgcheck's container was a
# REAL DEFECT the gates partly masked: a tmbstan built under
# StanHeaders >= 2.39 silently samples a standard normal instead of
# the model (dev/prior-dropping-investigation.md), which frm_sample()
# now refuses statically. The ungated correctness assertions are what
# caught it, which is why structural and exactness tests never take
# this switch. FRMTMB_SAMPLER_GATES=false remains in the pkgcheck
# workflow for reason (1) only.
sampler_gates_on <- function() {
  !identical(Sys.getenv("FRMTMB_SAMPLER_GATES"), "false")
}

# The skip that three of frmtmb's files each defined for themselves,
# defined once here now that they share a package.
skip_sampler <- function() {
  testthat::skip_if_not_installed("tmbstan")
  testthat::skip_if_not_installed("rstan")
}

# Copied from frmtmb's helper-reference.R rather than reached for: a
# three-line expectation is cheaper to duplicate than to make into API,
# which is the same call the package code made for `%||%` and the
# argument checks (see the (c) column of dev/draws-extraction.md).
expect_vector_equal <- function(x, y, tol) {
  expect_equal(length(x), length(y))
  expect_lt(max(abs(unname(x) - unname(y))), tol)
}

# The nonlinear loss-development model from frmtmb's brms-migration
# vignette, copied from its test-prior-compat.R with the blocks that
# use it. frmtmb keeps its own copy for the blocks that stayed there.
loss_data <- function(seed = 903) {
  set.seed(seed)
  AY <- factor(rep(1988:1997, each = 10))
  dev <- rep(seq(6, 114, by = 12), 10)
  ult_g <- 5000 + stats::rnorm(10, 0, 400)
  cum <- ult_g[as.integer(AY)] * (1 - exp(-(dev / 45)^1.3)) +
    stats::rnorm(100, 0, 120)
  data.frame(cum = cum, dev = dev, AY = AY)
}

loss_form <- function() {
  bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
     ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE) + gaussian()
}

# the vignette's own starting region; a nonlinear body this shaped has
# no useful default start
loss_start <- list(beta = c(5000, 1, 45))

# Data fixtures copied from the frmtmb test files these blocks came
# from. Each one stayed there too, because a block that did NOT move
# still uses it; copying is what keeps the two suites from depending on
# each other's private helpers.

# from tests/testthat/test-lkj.R
lkj_data <- function(seed = 77, n = 160L, ng = 16L) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n),
                   g = factor(rep(seq_len(ng), length.out = n)))
  u <- matrix(stats::rnorm(2 * ng), 2)
  dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + 0.8 * u[1, dd$g] +
                         0.4 * u[2, dd$g] * dd$x, 1)
  dd
}

# from tests/testthat/test-ordinal-fitted.R
ordfit_data <- function(seed, n = 250, tau = c(-0.8, 0.6), beta = 0.9,
                        levels = c("lo", "mid", "hi")) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n), z = stats::rnorm(n))
  eta <- beta * dd$x
  p <- cbind(stats::plogis(tau[1] - eta),
             stats::plogis(tau[2] - eta) - stats::plogis(tau[1] - eta),
             1 - stats::plogis(tau[2] - eta))
  dd$y <- factor(levels[apply(p, 1L, function(pr) sample(3L, 1L,
                                                         prob = pr))],
                 levels = levels, ordered = TRUE)
  dd
}

# from tests/testthat/test-priors-bounds-grcov.R
sim_lmm <- function(seed = 301, n = 150, ng = 15) {
  set.seed(seed)
  dd <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng),
                                                length.out = n)))
  dd$y <- rnorm(n, 1 + 0.5 * dd$x + rnorm(ng, 0, 0.7)[dd$g], 1)
  dd
}

# from tests/testthat/test-hmm.R
sim_hmm <- function(N, Tl, G, mu, sigma, seed) {
  set.seed(seed)
  K <- nrow(G)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl); s[1L] <- sample.int(K, 1L)
    for (t in seq_len(Tl)[-1L]) {
      s[t] <- sample.int(K, 1L, prob = G[s[t - 1L], ])
    }
    data.frame(id = id, t = seq_len(Tl), state = s,
               y = stats::rnorm(Tl, mu[s], sigma[s]))
  }))
}

# from tests/testthat/test-lca.R
# simulated K = 2, J = 4 binary-item data with well-separated classes
sim_lca_data <- function(seed = 5, n = 400, p1 = 0.35,
                         pr = rbind(c(0.90, 0.85, 0.20, 0.75),
                                    c(0.20, 0.15, 0.85, 0.25))) {
  set.seed(seed)
  cl <- rbinom(n, 1, p1) + 1
  J <- ncol(pr)
  Y <- matrix(0L, n, J)
  for (j in seq_len(J)) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
  dd <- data.frame(x = rnorm(n))
  dd$Y <- Y
  list(dd = dd, cl = cl)
}

# from tests/testthat/test-review-v29.R
v29_ordinal_data <- function(seed, n = 250, tau = c(-0.8, 0.6),
                             beta = 0.9) {
  set.seed(seed)
  dd <- data.frame(x = stats::rnorm(n))
  eta <- beta * dd$x
  p <- cbind(stats::plogis(tau[1] - eta),
             stats::plogis(tau[2] - eta) - stats::plogis(tau[1] - eta),
             1 - stats::plogis(tau[2] - eta))
  dd$y <- factor(apply(p, 1L, function(pr) sample(3L, 1L, prob = pr)),
                 levels = 1:3, ordered = TRUE)
  dd
}

# from tests/testthat/test-hmm.R
G2 <- matrix(c(0.9, 0.1, 0.25, 0.75), 2, 2, byrow = TRUE)

# from tests/testthat/test-setprior.R
fit_sp <- local({
  set.seed(401)
  dd <- data.frame(x = rnorm(150), f = factor(rep(c("a", "b"), 75)),
                   g = factor(rep(1:15, 10)))
  dd$y <- rnorm(150, 1 + 0.5 * dd$x + (dd$f == "b") +
                  rnorm(15, 0, 0.7)[dd$g], 1)
  frm(bf(y ~ x + f + (1 | g), sigma ~ x) + gaussian(), data = dd)
})

# Punch round 1, minor 8: extensions whose compat rows call autoscale
# "untested" now get it by default below sd 1e-3. For each, a covariate
# on one of the family's own dpars at scale 1 and at 1e-6: does the
# default engage, and does it reach the scale-1 logLik? autoscale =
# FALSE shown beside it.
#   PREDFIX_ARM=lane Rscript dev/predfix-p1-ext.R > dev/predfix-log/p1-ext.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
suppressPackageStartupMessages({
  library(frmtmb.latent); library(frmtmb.coupling); library(frmtmb.learn)
})
q <- function(expr) suppressWarnings(suppressMessages(expr))
ll <- function(f) if (inherits(f, "try-error")) NA else as.numeric(logLik(f))
run <- function(label, make_fit, d, col) {
  x0 <- d[[col]]
  ref <- ll(try(q(make_fit(d, NULL)), silent = TRUE))
  d[[col]] <- x0 * 1e-6
  fd <- try(q(make_fit(d, NULL)), silent = TRUE)
  ff <- try(q(make_fit(d, FALSE)), silent = TRUE)
  cat(sprintf("%-40s scale-1 %14.6f | 1e-6 default %14.6f (engaged %s, short %.3g) | FALSE short %.3g\n",
              label, ref, ll(fd),
              if (inherits(fd, "try-error")) "ERR" else !is.null(fd$par_units),
              ref - ll(fd), ref - ll(ff)))
  if (inherits(fd, "try-error")) cat("   error:", conditionMessage(attr(fd, "condition")), "\n")
}
ctl <- function(a) frmtmb_control(autoscale = a)

# lca: latent class regression, the covariate on the class logit
for (seed in 1:3) {
  set.seed(seed)
  n <- 400
  x <- rnorm(n)
  cl <- 1L + rbinom(n, 1, plogis(-0.8 + 1.2 * x))
  pr <- rbind(c(0.90, 0.85, 0.20, 0.75), c(0.15, 0.10, 0.90, 0.20))
  Y <- matrix(0L, n, 4)
  for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
  dd <- data.frame(x = x)
  dd$Y <- Y
  run(paste("lca(K = 2), Y ~ x, seed", seed),
      function(d, a) frm(bf(Y ~ x), family = lca(K = 2), data = d,
                         control = ctl(a)), dd, "x")
}

# cross_wishart: a covariate on the coherence dpar
for (seed in 1:3) {
  set.seed(seed)
  N <- 200
  n <- 8
  z <- rnorm(N)
  rows <- lapply(seq_len(N), function(i) {
    coh <- plogis(-0.2 + 0.6 * z[i])
    a <- sqrt(1.4); cm <- sqrt(coh * 0.8); b <- sqrt((1 - coh) * 0.8)
    z1 <- complex(real = rnorm(n, 0, sqrt(0.5)), imaginary = rnorm(n, 0, sqrt(0.5)))
    z2 <- complex(real = rnorm(n, 0, sqrt(0.5)), imaginary = rnorm(n, 0, sqrt(0.5)))
    d1 <- a * z1
    d2 <- cm * complex(modulus = 1, argument = -0.6) * z1 + b * z2
    cr <- sum(d1 * Conj(d2))
    c(sum(Mod(d1)^2), sum(Mod(d2)^2), Re(cr), Im(cr))
  })
  m <- do.call(rbind, rows)
  dd <- data.frame(w11 = m[, 1], w22 = m[, 2], w12r = m[, 3], w12i = m[, 4],
                   n = n, z = z)
  run(paste("cross_wishart, coh ~ z, seed", seed),
      function(d, a) frm(bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
                            pow2 ~ 1, coh ~ z, phase ~ 1),
                         family = cross_wishart(), data = d,
                         control = ctl(a)), dd, "z")
}

# frmtmb.learn: bandit2arm_delta, a subject covariate on tau
for (seed in 1:3) {
  d <- frm_task_design("bandit2arm", n_subject = 30, n_trial = 100,
                       seed = seed)
  set.seed(seed)
  zs <- rnorm(30)
  i <- as.integer(d$id)
  d$z <- zs[i]
  d$choice <- frm_task_simulate(
    bandit2arm_delta(subject = id, trial = trial), d,
    pars = list(alpha = 0.35, tau = exp(log(3) + 0.3 * zs[i])),
    seed = seed)[[1L]]$choice
  run(paste("bandit2arm_delta, tau ~ z, seed", seed),
      function(dd, a) frm(bf(choice | reward(pay1, pay2) ~ 1, tau ~ z),
                          family = bandit2arm_delta(subject = id,
                                                    trial = trial),
                          data = dd, control = ctl(a)), d, "z")
}

# hmm(K = 2, gaussian()): a covariate on the state means' shared mu
# predictor; the data carry a real slope so the column matters
for (seed in 1:3) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.25, 0.75), 2, 2, byrow = TRUE)
  dd <- do.call(rbind, lapply(1:25, function(id) {
    s <- integer(24); s[1L] <- sample.int(2, 1L)
    for (t in 2:24) s[t] <- sample.int(2, 1L, prob = G[s[t - 1L], ])
    x <- rnorm(24)
    data.frame(id = id, t = 1:24, x = x,
               y = rnorm(24, c(0, 3)[s] + 0.5 * x, 0.6))
  }))
  run(paste("hmm(K = 2, gaussian()), y ~ x, seed", seed),
      function(d, a) frm(bf(y ~ x), family = hmm(K = 2, gaussian(), time = t,
                                                 group = id),
                         data = d, control = ctl(a)), dd, "x")
}

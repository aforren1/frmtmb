# Attack 9: residuals(type = "osa") on a latent fit, which the lane checked only
# to run. On data from the model the residuals must be near iid N(0, 1). The
# reference is the exact one-step predictive by the forward recursion on a grid.
# Also unsorted rows.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
sim <- function(seed, fam, G = 15, K = 10, s = .6, ar = .6, b0 = 1) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = seq_len(G))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  e <- as.vector(replicate(G, arima.sim(list(ar = ar), K, sd = s)))
  eta <- b0 + .4 * d$x + e
  d$y <- if (fam == "poisson") rpois(nrow(d), exp(eta)) else rgamma(nrow(d), 5, 5 / exp(eta))
  d$g <- factor(d$g)
  d
}
exact_osa <- function(d, fam, b0, b1, s, ar, shape = NA, ng = 161) {
  th <- ar / sqrt(1 - ar^2)
  msd <- s * sqrt(1 + th^2)
  u <- seq(-6.5, 6.5, length.out = ng)
  h <- u[2] - u[1]
  grid <- u * msd
  Tm <- outer(grid, grid, function(a, b) dnorm(b, ar * a, s)) * (h * msd)
  out <- numeric(nrow(d))
  cdf <- function(y, mu) if (fam == "poisson") (ppois(y - 1, mu) + ppois(y, mu)) / 2 else pgamma(y, shape, shape / mu)
  dens <- function(y, mu) if (fam == "poisson") dpois(y, mu) else dgamma(y, shape, shape / mu)
  for (idx in split(seq_len(nrow(d)), d$g)) {
    idx <- idx[order(d$t[idx])]
    f <- dnorm(u) * h
    for (k in seq_along(idx)) {
      i <- idx[k]
      if (k > 1L) f <- as.vector(crossprod(f, Tm))
      f <- f / sum(f)
      mu <- exp(b0 + b1 * d$x[i] + grid)
      out[i] <- qnorm(pmin(pmax(sum(f * cdf(d$y[i], mu)), 1e-12), 1 - 1e-12))
      f <- f * dens(d$y[i], mu)
    }
  }
  out
}
acf1 <- function(r, d) {
  o <- order(d$g, d$t)
  r <- r[o]; g <- d$g[o]
  keep <- g[-1] == g[-length(g)]
  cor(r[-1][keep], r[-length(r)][keep], use = "complete.obs")
}
rows <- list()
for (fam in c("poisson", "Gamma")) {
  for (seed in 1:20) {
    d <- sim(seed, fam)
    fm <- if (fam == "poisson") poisson() else Gamma(link = "log")
    f <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = fm, data = d))
    th <- f$estimates$theta
    shape <- if (fam == "Gamma") exp(f$opt$par[names(f$opt$par) == "betad"]) else NA
    ar <- th[2] / sqrt(1 + th[2]^2)
    if (abs(ar) > .98) next
    t0 <- proc.time()[[3]]
    ro <- msg_of(residuals(f, type = "osa"))
    tosa <- proc.time()[[3]] - t0
    if (!is.null(ro$err)) {
      cat("ERR", seed, fam, ro$err, "\n")
      next
    }
    r <- as.matrix(ro$val)[, 1]
    x <- exact_osa(d, fam, f$estimates$beta[1], f$estimates$beta[2], exp(th[1]), ar, shape)
    rows[[length(rows) + 1]] <- data.frame(fam = fam, seed = seed, tosa = tosa, mean_frm = mean(r), sd_frm = sd(r), acf_frm = acf1(r, d),
      mean_ex = mean(x), sd_ex = sd(x), acf_ex = acf1(x, d), cor_frm_ex = cor(r, x),
      maxabs_diff = max(abs(r - x)), ks_frm = suppressWarnings(ks.test(r, "pnorm")$p.value))
    # unsorted rows
    if (seed <= 5) {
      set.seed(seed + 1000)
      ds <- d[sample(nrow(d)), ]
      fs <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = fm, data = ds))
      rs <- as.matrix(residuals(fs, type = "osa"))[, 1]
      rows[[length(rows)]]$unsorted_sd <- sd(rs)
      rows[[length(rows)]]$unsorted_mean <- mean(rs)
    } else {
      rows[[length(rows)]]$unsorted_sd <- NA
      rows[[length(rows)]]$unsorted_mean <- NA
    }
  }
}
o <- do.call(rbind, rows)
print(o, digits = 3, row.names = FALSE)
print(aggregate(cbind(mean_frm, sd_frm, acf_frm, mean_ex, sd_ex, acf_ex, cor_frm_ex, maxabs_diff, ks_frm, tosa) ~ fam, o, mean), digits = 3)
print(residuals(f, type = "osa")[1:2, ])

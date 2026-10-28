# Reviewer probe 2a: is frmtmb's prior = "flat" density on a correlation
# block Stan's implicit uniform on cholesky_factor_corr, for d = 2, 3, 4,
# on us() blocks, rescor and unstr()?
#  (1) importance check: draw t row by row from the exact density that
#      the derivation says flat-on-L induces on frmtmb's chart (a
#      bivariate-t-like (1 + |t|^2)^(-(n + 2) / 2) per row, multivariate
#      t with 2 df and scale I / 2); the log ratio of frmtmb's RESOLVED
#      density to it must be constant. lkj(1) for contrast.
#  (2) independent reference: Stan's draws of `cholesky_factor_corr[D] L`
#      with no model block (brms's program with "" on cor), compared on
#      E[r_ij^2] and P(r_ij > 0.5) with the frmtmb draws of (1).
# Usage: Rscript dev/review-reslot-p2a.R after. Seeds 1 (data), 21
# (draws), 31 (Stan).
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
suppressPackageStartupMessages(library(rstan))
set.seed(1)
G <- 8L; n <- 12L
dd <- data.frame(g = factor(rep(seq_len(G), each = n)), t = rep(seq_len(n), G))
dd$x <- rnorm(nrow(dd)); dd$z <- rnorm(nrow(dd)); dd$w <- rnorm(nrow(dd))
dd$y <- 1 + 0.5 * dd$x + rnorm(G, 0, 0.7)[dd$g] + rnorm(nrow(dd))
dd$y2 <- 0.5 * dd$y + rnorm(nrow(dd)); dd$y3 <- 0.3 * dd$y + rnorm(nrow(dd))
dd$y4 <- rnorm(nrow(dd))
dd$ts <- (dd$t - 1) %% 4 + 1
dd$gs <- factor(paste(dd$g, (dd$t - 1) %/% 4))

lkj_entry <- function(fit, prior) {
  r <- suppressWarnings(suppressMessages(
    frmtmb.sample:::sample_resolve_priors(fit, prior)))
  e <- Filter(function(e) identical(e$dist$kind, "lkj"), r$ri$entries)
  stopifnot(length(e) == 1L)
  e[[1L]]
}
# frmtmb's resolved log density at internal t (the entry alone)
ld_entry <- function(e, t) {
  frmtmb:::prior_logdens(t, e$dist, e$scale)
}
rows_of <- function(d) {
  ii <- row(diag(d))[lower.tri(diag(d))]
  lapply(2:d, function(i) which(ii == i))
}
# exact draw from the claimed induced density, and its log density
draw_t <- function(d) {
  t <- numeric(d * (d - 1) / 2); lq <- 0
  for (r in rows_of(d)) {
    k <- length(r)
    z <- rnorm(k) / sqrt(2)
    s <- sqrt(rchisq(1, 2) / 2)
    v <- z / s
    t[r] <- v
    lq <- lq + lgamma((2 + k) / 2) - lgamma(1) - k / 2 * log(2 * pi) +
      k / 2 * log(2) - (2 + k) / 2 * log(1 + sum(v^2))
  }
  list(t = t, lq = lq)
}
cor_of <- function(t, d) {
  C <- frmtmb:::us_chol_cor(t, d)
  C[lower.tri(C)]
}
stan_code <- "data { int<lower=2> D; } parameters { cholesky_factor_corr[D] L; }
generated quantities { matrix[D, D] R = multiply_lower_tri_self_transpose(L); }"
sm <- rstan::stan_model(model_code = stan_code)
stan_draws <- function(d) {
  f <- rstan::sampling(sm, data = list(D = d), chains = 4, iter = 6000,
                       warmup = 1000, seed = 31, refresh = 0)
  R <- rstan::extract(f, "R")$R
  matrix(t(apply(R, 1, function(m) m[lower.tri(m)])), nrow = dim(R)[1])
}

cases <- list(
  list(label = "us d=2 (1 + x | g)", f = bf(y ~ (1 + x | g)), d = 2),
  list(label = "us d=3 (1 + x + z | g)", f = bf(y ~ (1 + x + z | g)), d = 3),
  list(label = "us d=4 (1 + x + z + w | g)", f = bf(y ~ (1 + x + z + w | g)),
       d = 4),
  list(label = "rescor d=3", f = bf(y ~ x) + bf(y2 ~ x) + bf(y3 ~ x) +
         set_rescor(TRUE), d = 3),
  list(label = "rescor d=4", f = bf(y ~ x) + bf(y2 ~ x) + bf(y3 ~ x) +
         bf(y4 ~ x) + set_rescor(TRUE), d = 4),
  list(label = "unstr d=4", f = bf(y ~ x + unstr(ts, gs)), d = 4))
set.seed(21)
N <- 20000L
for (cs in cases) {
  fit <- suppressWarnings(suppressMessages(
    frm(cs$f, data = dd, dry_run = "objective")))
  ef <- lkj_entry(fit, "flat")
  e1 <- lkj_entry(fit, NULL)          # the default lkj(1)
  d <- cs$d
  lw_f <- lw_1 <- numeric(200)
  for (j in 1:200) {
    dr <- draw_t(d)
    lw_f[j] <- ld_entry(ef, dr$t) - dr$lq
    lw_1[j] <- ld_entry(e1, dr$t) - dr$lq
  }
  cat(sprintf("== %s: idx %s; log weight range: flat %.3g, lkj(1) %.3g\n",
              cs$label, paste(ef$idx, collapse = ","), diff(range(lw_f)),
              diff(range(lw_1))))
}
cat("\nMoments of r_ij (lower triangle, column-major): frmtmb flat draws",
    "(exact, N =", N, ") vs Stan cholesky_factor_corr with no model\n")
for (d in 2:4) {
  ft <- t(replicate(N, cor_of(draw_t(d)$t, d)))
  if (d == 2) ft <- t(ft)
  st <- stan_draws(d)
  # lkj(1) reference: rejection-free via onion is overkill; sample the
  # Stan program with lkj_corr_cholesky(1) is not needed: E[r^2] under
  # lkj(1) is 1 / (d + 1) for every entry
  for (k in seq_len(ncol(st))) {
    se <- sqrt(var(st[, k]^2) / nrow(st) * 4 + var(ft[, k]^2) / N)
    cat(sprintf("  d=%d r%d: E[r^2] frmtmb %.4f Stan %.4f (diff/se %.2f) lkj(1) %.4f | P(r>0.5) frmtmb %.4f Stan %.4f\n",
                d, k, mean(ft[, k]^2), mean(st[, k]^2),
                (mean(ft[, k]^2) - mean(st[, k]^2)) / se, 1 / (d + 1),
                mean(ft[, k] > 0.5), mean(st[, k] > 0.5)))
  }
}

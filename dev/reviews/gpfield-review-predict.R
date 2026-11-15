# Review of lane gpfield, item 1919: is predict()'s draw at a new position
# the exact Gaussian conditional given the replicate's own field?
# Per replicate s we record fs (the replicate's parameters and field) and
# the mu predict() computed at the new positions, then compute OUR OWN
# conditional mean and covariance of f(x*) given f_s at the training
# positions at theta_s (kernel written here, not frmtmb's), and whiten.
# Usage: Rscript gpfield-review-predict.R <lane|base> <cov> <ndraws>
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]; cv <- a[2]; nd_ <- as.integer(a[3])
gf_libs(arm)
ns <- asNamespace("frmtmb")
kern <- function(r, ls, cv) {
  s <- r / ls
  switch(cv, exp_quad = exp(-s^2 / 2), exponential = exp(-s),
         matern32 = (1 + sqrt(3) * s) * exp(-sqrt(3) * s),
         matern52 = (1 + sqrt(5) * s + 5 / 3 * s^2) * exp(-sqrt(5) * s))
}
set.seed(9200)
n <- 30
d <- data.frame(x = sort(stats::runif(n, 0, 10)))
d$y <- 2 * sin(d$x) + stats::rnorm(n, 0, 0.2)
fml <- if (cv == "exp_quad") y ~ gp(x) else
  as.formula(sprintf('y ~ gp(x, cov = "%s")', cv))
q <- gf_fit(bf(fml), data = d, family = gaussian())
fit <- q$fit
cat("fit conv", fit$opt$convergence, "warn", length(q$warn), q$warn, "\n")
bk <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)[[1]]
lp <- fit$frame$linpreds[[ns$gp_block_lp(bk, fit$frame)]]
xs <- c(5.5, 10.5, 12, 12.3, 15, 60)
newd <- data.frame(x = xs)
rec <- new.env()
rec$fs <- list(); rec$mu <- list()
orig <- ns$predict_dpar_values
gf_put("predict_dpar_values", function(fs, ...) {
  out <- orig(fs, ...)
  k <- length(rec$fs) + 1L
  rec$fs[[k]] <- fs$estimates
  rec$mu[[k]] <- out
  out
})
set.seed(if (length(a) > 3) as.integer(a[4]) else 1L)
t0 <- proc.time()[["elapsed"]]
p <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd_)
cat("predict wall", proc.time()[["elapsed"]] - t0, "ndraws", nd_, "\n")
gf_put("predict_dpar_values", orig)
cat("recorded", length(rec$fs), "\n")
str(rec$mu[[1]], max.level = 2)
xt <- d$x
Zg <- as.matrix(ns$lp_Z(fit, lp))[, bk$c_idx, drop = FALSE]
Z <- matrix(NA_real_, length(rec$fs), length(xs))
nsd <- numeric(0)
pr <- matrix(NA_real_, length(rec$fs), 2)
for (s in seq_along(rec$fs)) {
  est <- rec$fs[[s]]
  th <- est$theta[bk$theta_idx]
  sd2 <- exp(2 * th[1]); ls <- exp(th[2])
  cvec <- ns$expand_b(fit$frame, est$b, est$theta)
  ft <- as.vector(Zg %*% cvec[bk$c_idx])
  b0 <- est$beta[1]
  mu <- rec$mu[[s]]$mu
  if (is.null(mu)) mu <- rec$mu[[s]][[1]]
  mu <- as.numeric(mu)
  Ktt <- sd2 * (kern(abs(outer(xt, xt, "-")), ls, cv) + diag(1e-12, n))
  Kst <- sd2 * kern(abs(outer(xs, xt, "-")), ls, cv)
  Kss <- sd2 * kern(abs(outer(xs, xs, "-")), ls, cv)
  R <- chol(Ktt)
  # stable form: P = Kst R^-1, C = Kss - P P', m = P R^-T f
  P <- t(backsolve(R, t(Kst), transpose = TRUE))
  m <- as.vector(P %*% backsolve(R, ft, transpose = TRUE))
  C <- Kss - tcrossprod(P)
  e <- mu - b0 - m
  Z[s, ] <- e / sqrt(pmax(diag(C), 1e-300))
  # correlation of the residual pair (12, 12.3)
  pr[s, ] <- c(C[3, 4] / sqrt(C[3, 3] * C[4, 4]), NA)
}
for (j in seq_along(xs)) {
  z <- Z[, j]
  cat(sprintf(paste0("COND %s x %5.1f | z mean %+.3f (se %.3f) sd %.3f | ",
                     "KS p %.3f | ndraw %d\n"),
              cv, xs[j], mean(z), sd(z) / sqrt(length(z)), sd(z),
              suppressWarnings(stats::ks.test(z, "pnorm")$p.value),
              length(z)))
}
cat(sprintf("COND %s pair (12, 12.3) z-corr %.3f theory mean %.3f\n", cv,
            cor(Z[, 3], Z[, 4]), mean(pr[, 1])))
cat("SD of predict at xs:", round(apply(p, 2, sd), 3), "\n")

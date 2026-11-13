# Review of lane gpstop: the sweep's student rows where the relinearization
# moved the fit (logis_a n40 seed 2, pow_b n40 seed 2), both arms: the
# inner Hessian's smallest eigenvalue, the Laplace value, and the marginal
# by importance sampling (t df 4, floors 1 and 0.1, seeds 9001/9002).
# Data: gpstop-review-sweep2.R mk(),
# set.seed(81100 + 1000*4 + 100*seed + n + 10*body index), y = eta + 0.1 t_3.
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
mk <- function(body, n, seed) {
  set.seed(81100 + 1000 * 4 + 100 * seed + n + 10 * match(body, bodies))
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  w <- round(runif(n, 0.2, 5), 3)
  f <- if (seed %% 2) sin(x) else 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  eta <- switch(body,
    logis_a = (3 + fs) / (1 + exp(-1.5 * z)),
    pow_b = w^(0.5 + 0.3 * fs))
  data.frame(x, z, w, y = eta + 0.1 * rt(n, 3))
}
is_ll <- function(obj, par, floor, seed, ndraw = 20000L, df = 4) {
  v <- -obj$fn(par)
  env <- obj$env
  full <- env$last.par
  r <- env$random
  H <- as.matrix(env$spHess(full, random = TRUE))
  H <- 0.5 * (H + t(H))
  e <- eigen(H, symmetric = TRUE)
  P <- e$vectors %*% diag(pmax(e$values, floor)) %*% t(e$vectors)
  R <- chol(0.5 * (P + t(P)))
  k <- length(r)
  set.seed(seed)
  z <- matrix(rnorm(ndraw * k), ndraw, k)
  s <- sqrt(rchisq(ndraw, df) / df)
  dev <- t(backsolve(R, t(z))) / s
  lq <- lgamma((df + k) / 2) - lgamma(df / 2) - k / 2 * log(df * pi) +
    sum(log(diag(R))) - (df + k) / 2 * log1p(rowSums(z^2) / s^2 / df)
  u0 <- full[r]
  lj <- vapply(seq_len(ndraw), function(i) {
    x <- full
    x[r] <- u0 + dev[i, ]
    -env$f(x, order = 0)
  }, 0)
  lw <- lj - lq
  m <- max(lw)
  w <- exp(lw - m)
  c(lap = v, is = m + log(mean(w)), ess = sum(w)^2 / sum(w^2),
    ev_min = min(e$values))
}
specs <- list(
  logis_a = bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE),
  pow_b = bf(y ~ a * w^b, a ~ 1, b ~ gp(x), nl = TRUE))
for (body in names(specs)) {
  d <- mk(body, 40L, 2L)
  for (arm in c("lane", "norelin")) {
    q <- rv_run(paste0(body, "_student_n40_s2"), arm, specs[[body]], data = d,
                family = student())
    f <- q$fit
    for (cfg in list(c(1, 9001), c(1, 9002), c(0.1, 9001))) {
      v <- is_ll(f$obj, f$opt$par, cfg[1], cfg[2])
      cat(sprintf(paste0("DEG2IS %s %s floor %g seed %d laplace %.4f is %.4f ",
                         "ess %.0f inner ev_min %.3g\n"),
                  body, arm, cfg[1], cfg[2], v[["lap"]], v[["is"]], v[["ess"]],
                  v[["ev_min"]]))
    }
    flush.console()
  }
}
cat("DONE\n")

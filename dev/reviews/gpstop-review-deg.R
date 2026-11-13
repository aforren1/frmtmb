# Review of lane gpstop: the sweep's logis_a bernoulli n40 seed 2 row, where
# the lane keeps log length scale -323 and 0.81.0's route log sd -40.6.
# Both arms; the limit model with an observation-level iid effect in `a`;
# the marginal by importance sampling at each optimum (t df 4, eigenvalues
# floored at 1 and at 0.1, set.seed 9001/9002); the warning text of the
# mm_K bernoulli n120 seed 1 row.
# Data: gpstop-review-sweep.R mk(), set.seed(81100 + 1000*3 + 100*2 + 40 + 30)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
mk <- function(body, fam, n, seed) {
  set.seed(81100 + 1000 * match(fam, c("gaussian", "poisson", "bernoulli")) +
             100 * seed + n + 10 * match(body, bodies))
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  w <- round(runif(n, 0.2, 5), 3)
  f <- if (seed %% 2) sin(x) else 2 * sin(3 * x) + x
  fs <- f / max(abs(f))
  eta <- switch(body,
    logis_a = (3 + fs) / (1 + exp(-1.5 * z)),
    mm_K = 3 * w / (1 + 0.6 * fs + w))
  y <- rbinom(n, 1, plogis(1.5 * (eta - median(eta)) / sd(eta)))
  data.frame(x, z, w, y)
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
d <- mk("logis_a", "bernoulli", 40L, 2L)
fml <- bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE)
for (arm in c("lane", "norelin")) {
  q <- rv_run("logis_bern40s2", arm, fml, data = d, family = bernoulli())
  f <- q$fit
  cat("DEGP", arm, paste(names(f$opt$par), signif(f$opt$par, 5),
                         collapse = " "), "\n")
  for (cfg in list(c(1, 9001), c(1, 9002), c(0.1, 9001))) {
    v <- is_ll(f$obj, f$opt$par, cfg[1], cfg[2])
    cat(sprintf("DEGIS %s floor %g seed %d laplace %.4f is %.4f ess %.0f ev_min %.3g\n",
                arm, cfg[1], cfg[2], v[["lap"]], v[["is"]], v[["ess"]],
                v[["ev_min"]]))
  }
  s <- tryCatch(capture.output(print(summary(f)))[1:40],
                error = function(e) conditionMessage(e))
  cat(paste("DEGS", arm, s), sep = "\n")
  flush.console()
}
rv_arm("lane")
d$obs <- factor(seq_len(nrow(d)))
q <- rv_fit(bf(y ~ a / (1 + exp(-b * z)), a ~ 1 + (1 | obs), b ~ 1,
               nl = TRUE), data = d, family = bernoulli())
cat(sprintf("DEGLIM iid obs effect ll %.4f conv %d nwarn %d | %s\n",
            as.numeric(logLik(q$fit)), q$fit$opt$convergence, length(q$warn),
            paste(names(q$fit$opt$par), signif(q$fit$opt$par, 5),
                  collapse = " ")))
q <- rv_fit(bf(y ~ a / (1 + exp(-b * z)), a ~ 1, b ~ 1, nl = TRUE),
            data = d, family = bernoulli())
cat(sprintf("DEGLIM no field ll %.4f conv %d nwarn %d\n",
            as.numeric(logLik(q$fit)), q$fit$opt$convergence, length(q$warn)))
# the extra warning on mm_K bernoulli n120 seed 1
d2 <- mk("mm_K", "bernoulli", 120L, 1L)
q <- rv_fit(bf(y ~ V * w / (K + w), V ~ 1, K ~ gp(x), nl = TRUE), data = d2,
            family = bernoulli())
cat("DEGW mmK120s1 lane", q$warn, sep = "\n  ")
cat("DONE\n")

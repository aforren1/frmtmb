# Lane gpstop, item 2: is the reference's higher student() maximum a
# maximum of the model, or of the Laplace approximation alone? At the
# default's and the reference's outer optima: the inner Hessian's smallest
# eigenvalues, and the marginal log-likelihood by importance sampling with a
# multivariate t (df 4) proposal at the inner mode and inverse Hessian,
# 20000 draws, set.seed(9001). Args: arm, rows as case:curve:n:seed:rho0:sd0.
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
gs_libs(a[1])
rows <- strsplit(a[2], ",")[[1]]
ndraw <- if (length(a) > 2) as.integer(a[3]) else 20000L
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x,
               wig = function(x) sin(7 * x) + 0.2 * x,
               two = function(x) sin(x) + 0.4 * sin(10 * x),
               smooth = function(x) sin(x), step = function(x) 2 * (x > 5) - 1)
is_ll <- function(obj, par, df = 4) {
  env <- obj$env
  lap <- -obj$fn(par)
  full <- env$last.par
  r <- env$random
  H <- as.matrix(env$spHess(full, random = TRUE))
  H <- 0.5 * (H + t(H))
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  u0 <- full[r]
  R <- chol(H)
  k <- length(r)
  set.seed(9001)
  z <- matrix(rnorm(ndraw * k), ndraw, k)
  s <- sqrt(rchisq(ndraw, df) / df)
  # u = u0 + R^-1 z / s, a t draw with scale H^-1
  dev <- t(backsolve(R, t(z))) / s
  lq <- lgamma((df + k) / 2) - lgamma(df / 2) - k / 2 * log(df * pi) +
    sum(log(diag(R))) - (df + k) / 2 * log1p(rowSums(z^2) / s^2 / df)
  lj <- vapply(seq_len(ndraw), function(i) {
    x <- full
    x[r] <- u0 + dev[i, ]
    -env$f(x, order = 0)
  }, 0)
  lw <- lj - lq
  m <- max(lw)
  w <- exp(lw - m)
  ess <- sum(w)^2 / sum(w^2)
  c(laplace = lap, is = m + log(mean(w)), ess = ess, ev_min = min(ev),
    ev_2 = sort(ev)[2], ev_max = max(ev))
}
for (rw in rows) {
  p <- strsplit(rw, ":")[[1]]
  cs <- p[1]; cn <- p[2]; n <- as.integer(p[3]); seed <- as.integer(p[4])
  rho0 <- as.numeric(p[5]); sd0 <- as.numeric(p[6])
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  d <- data.frame(x = x, y = curves[[cn]](x) + 0.02 * rt(n, 5))
  fd <- fe_fit(y ~ gp(x), data = d, family = student())$fit
  fr <- fe_fit(y ~ gp(x), data = d, family = student(),
               start = list(theta = c(log(sd0), log(rho0))))$fit
  for (wh in c("default", "ref")) {
    f <- if (wh == "default") fd else fr
    v <- is_ll(f$obj, f$opt$par)
    cat(sprintf(paste0("STIS %s %s %d %d %-7s laplace %.4f is %.4f ess %.0f ",
                       "of %d | inner H eigen min %.3g second %.3g max %.3g\n"),
                cs, cn, n, seed, wh, v[["laplace"]], v[["is"]], v[["ess"]],
                ndraw, v[["ev_min"]], v[["ev_2"]], v[["ev_max"]]))
    flush.console()
  }
}
cat("DONE\n")

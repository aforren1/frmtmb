# Lane gpstop, item 2: dev/gpstop-stis.R with a proposal that stays
# usable at a near-singular inner Hessian. The field's coefficients have a
# unit prior, so the proposal's precision takes each eigenvalue of the
# inner Hessian floored at 1 (a no-op where the Hessian holds the prior's
# curvature); multivariate t, df 4, set.seed(9001). Also the share of the
# Laplace value that the smallest eigenvalue carries, -log(ev_min) / 2.
# Args: arm, rows as case:curve:n:seed:rho0:sd0, draws.
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
  e <- eigen(H, symmetric = TRUE)
  ev <- e$values
  P <- e$vectors %*% diag(pmax(ev, 1)) %*% t(e$vectors)
  u0 <- full[r]
  R <- chol(0.5 * (P + t(P)))
  k <- length(r)
  set.seed(9001)
  z <- matrix(rnorm(ndraw * k), ndraw, k)
  s <- sqrt(rchisq(ndraw, df) / df)
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
  # the estimate on each half of the draws, for its spread
  h <- seq_len(ndraw) <= ndraw / 2
  half <- c(m + log(mean(w[h])), m + log(mean(w[!h])))
  c(laplace = lap, is = m + log(mean(w)), h1 = half[1], h2 = half[2],
    ess = sum(w)^2 / sum(w^2), ev_min = min(ev), share = -log(min(ev)) / 2)
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
    cat(sprintf(paste0("STIS2 %s %s %d %d %-7s laplace %.4f is %.4f ",
                       "(halves %.4f %.4f) ess %.0f of %d | ev_min %.3g, ",
                       "-log(ev_min)/2 %.3f\n"),
                cs, cn, n, seed, wh, v[["laplace"]], v[["is"]], v[["h1"]],
                v[["h2"]], v[["ess"]], ndraw, v[["ev_min"]], v[["share"]]))
    flush.console()
  }
}
cat("DONE\n")

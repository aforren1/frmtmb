# Review of lane gpstop, item 5033, st80 (student step 80, seed 1): fits
# from points on the segment between the default optimum and the
# reference, every outer component started there. Is there a regular
# maximum (code 0, no warning, healthy inner Hessian) above the default?
# Args: n, rho0:sd0 of the reference, t values (comma), draws.
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
n <- as.integer(a[1])
rs <- as.numeric(strsplit(a[2], ":")[[1]])
ts <- as.numeric(strsplit(a[3], ",")[[1]])
ndraw <- as.integer(a[4])
set.seed(74000 + 1 * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
d <- data.frame(x = x, y = 2 * (x > 5) - 1 + 0.02 * rt(n, 5))
fd <- rv_fit(y ~ gp(x), data = d, family = student())$fit
fr <- rv_fit(y ~ gp(x), data = d, family = student(),
             start = list(theta = c(log(rs[2]), log(rs[1]))))$fit
p0 <- fd$opt$par
p1 <- fr$opt$par
pn <- names(p0)
inner <- function(obj, par) {
  v <- -obj$fn(par)
  H <- as.matrix(obj$env$spHess(obj$env$last.par, random = TRUE))
  list(lap = v, full = obj$env$last.par, H = 0.5 * (H + t(H)))
}
is_ll <- function(obj, st, floor = 1, df = 4, seed = 9001) {
  env <- obj$env
  r <- env$random
  e <- eigen(st$H, symmetric = TRUE)
  P <- e$vectors %*% diag(pmax(e$values, floor)) %*% t(e$vectors)
  R <- chol(0.5 * (P + t(P)))
  k <- length(r)
  set.seed(seed)
  z <- matrix(rnorm(ndraw * k), ndraw, k)
  s <- sqrt(rchisq(ndraw, df) / df)
  dev <- t(backsolve(R, t(z))) / s
  lq <- lgamma((df + k) / 2) - lgamma(df / 2) - k / 2 * log(df * pi) +
    sum(log(diag(R))) - (df + k) / 2 * log1p(rowSums(z^2) / s^2 / df)
  u0 <- st$full[r]
  xx <- st$full
  lj <- vapply(seq_len(ndraw), function(i) {
    xx[r] <- u0 + dev[i, ]
    -env$f(xx, order = 0)
  }, 0)
  lw <- lj - lq
  m <- max(lw)
  w <- exp(lw - m)
  c(is = m + log(mean(w)), ess = sum(w)^2 / sum(w^2))
}
report <- function(tag, q) {
  f <- q$fit
  st <- inner(f$obj, f$opt$par)
  ev <- min(eigen(st$H, symmetric = TRUE, only.values = TRUE)$values)
  g <- max(abs(f$obj$gr(f$opt$par)))
  st <- inner(f$obj, f$opt$par)
  v1 <- is_ll(f$obj, st, 1, 4, 9001)
  v2 <- is_ll(f$obj, st, 0.1, 3, 9002)
  cat(sprintf(paste0("ST2 n %d %s ll %.4f conv %d nwarn %d gmax %.2e ",
                     "ev_min %.3g | IS %.4f (ess %.0f) %.4f (ess %.0f) | ",
                     "theta %s | %s\n"),
              n, tag, -f$opt$objective, f$opt$convergence, length(q$warn), g,
              ev, v1[["is"]], v1[["ess"]], v2[["is"]], v2[["ess"]],
              paste(sprintf("%.4f", f$opt$par[pn == "theta"]), collapse = ","),
              paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
  flush.console()
}
report("default", list(fit = fd, warn = character(0)))
for (tt in ts) {
  pt <- p0 + tt * (p1 - p0)
  st <- list(beta = unname(pt[pn == "beta"]), betad = unname(pt[pn == "betad"]),
             theta = unname(pt[pn == "theta"]))
  q <- rv_fit(y ~ gp(x), data = d, family = student(), start = st)
  report(sprintf("from t %.2f", tt), q)
}
cat("DONE\n")

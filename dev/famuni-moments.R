# The post-fit moments of the three families (mean_fn, var_fn and
# trunc_mean_fn) against references computed another way: quadrature of
# the density, a long direct sum, and brms 2.23.0's own epred functions.
# Output: dev/famuni-out/moments.txt
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-300)
post <- function(fam) unclass(fam)[["post"]]

cat("\n=== frechet\n")
fam <- frechet()
P <- post(fam)
grid <- expand.grid(mu = c(0.5, 3, 40), nu = c(1.3, 2.5, 6, 30))
dens <- function(y, mu, nu) {
  s <- mu / gamma(1 - 1 / nu)
  exp(log(nu) - log(s) - (1 + nu) * log(y / s) - (y / s)^(-nu))
}
# split at 1 where the range reaches 0 or Inf, where the integrands of
# the GEV moments are steep; a tolerance of 1e-10 is what integrate()
# reaches on every piece here
q_int <- function(f, lo, hi) {
  one <- function(a, b) {
    # integrate() reports "roundoff error" on a few GEV pieces at 1e-10;
    # those are retried at 1e-8, still far inside what is being checked
    tryCatch(stats::integrate(f, a, b, rel.tol = 1e-10, abs.tol = 1e-13,
                              subdivisions = 5000L)$value,
             error = function(e) {
               stats::integrate(f, a, b, rel.tol = 1e-8,
                                subdivisions = 5000L)$value
             })
  }
  if (lo < 1 && hi > 1) one(lo, 1) + one(1, hi) else one(lo, hi)
}
m_ref <- mapply(function(mu, nu) {
  q_int(function(y) y * dens(y, mu, nu), 0, Inf)
}, grid$mu, grid$nu)
cat("mean_fn vs quadrature: max rel err",
    format(max(rel(P$mean_fn(as.list(grid), list()), m_ref)), digits = 3),
    "over", nrow(grid), "points\n")
ok <- grid$nu > 2
v_ref <- mapply(function(mu, nu) {
  q_int(function(y) (y - mu)^2 * dens(y, mu, nu), 0, Inf)
}, grid$mu[ok], grid$nu[ok])
cat("var_fn vs quadrature (nu > 2): max rel err",
    format(max(rel(P$var_fn(as.list(grid[ok, ]), list()), v_ref)),
           digits = 3), "over", sum(ok), "points; at nu <= 2:",
    paste(unique(P$var_fn(as.list(grid[!ok, ]), list())), collapse = " "),
    "\n")
lb <- 0.4 * grid$mu
ub <- 2 * grid$mu
tm_ref <- mapply(function(mu, nu, a, b) {
  q_int(function(y) y * dens(y, mu, nu), a, b) /
    q_int(function(y) dens(y, mu, nu), a, b)
}, grid$mu, grid$nu, lb, ub)
cat("trunc_mean_fn on (0.4 mu, 2 mu) vs quadrature: max rel err",
    format(max(rel(P$trunc_mean_fn(as.list(grid), list(), lb, ub), tm_ref)),
           digits = 3), "\n")
tm_ref1 <- mapply(function(mu, nu, a) {
  q_int(function(y) y * dens(y, mu, nu), a, Inf) /
    q_int(function(y) dens(y, mu, nu), a, Inf)
}, grid$mu, grid$nu, lb)
cat("trunc_mean_fn on (0.4 mu, Inf) vs quadrature: max rel err",
    format(max(rel(P$trunc_mean_fn(as.list(grid), list(), lb,
                                   rep(Inf, nrow(grid))), tm_ref1)),
           digits = 3), "\n")

cat("\n=== gen_extreme_value\n")
fam <- gen_extreme_value()
P <- post(fam)
grid <- expand.grid(mu = c(-1, 2), sigma = c(0.5, 2),
                    xi = c(-0.6, -0.2, -1e-8, 0, 1e-8, 0.1, 0.4, 0.8))
qf <- function(p, m, s, k) {
  e <- -log(p)
  if (k == 0) m - s * log(e) else m + s * (e^(-k) - 1) / k
}
# Y = Q(exp(-e)) for e ~ Exp(1): the moments as integrals over e, whose
# integrands are smoother at the ends than the quantile function on (0, 1)
qe <- function(e, m, s, k) {
  if (k == 0) m - s * log(e) else m + s * (e^(-k) - 1) / k
}
# and over t = log(e), on which the integrand is smooth and decays at
# both ends for xi < 1 (beyond t = -600 and 6 it is below 1e-40)
tw <- function(t) exp(t - exp(t))
m_ref <- mapply(function(m, s, k) {
  q_int(function(t) qe(exp(t), m, s, k) * tw(t), -600, 6)
}, grid$mu, grid$sigma, grid$xi)
mf <- P$mean_fn(as.list(grid), list())
cat("mean_fn vs quadrature over log(-log F): max rel err",
    format(max(rel(mf, m_ref)), digits = 3), "over", nrow(grid), "points\n")
ok <- grid$xi < 0.5
v_ref <- mapply(function(m, s, k, mm) {
  q_int(function(t) (qe(exp(t), m, s, k) - mm)^2 * tw(t), -600, 6)
}, grid$mu[ok], grid$sigma[ok], grid$xi[ok], m_ref[ok])
cat("var_fn vs quadrature (xi < 1/2): max rel err",
    format(max(rel(P$var_fn(as.list(grid[ok, ]), list()), v_ref)),
           digits = 3), "over", sum(ok), "points\n")
cat("mean at xi = 1 and 1.5:",
    format(P$mean_fn(list(mu = c(0, 0), sigma = c(1, 1), xi = c(1, 1.5)),
                     list())), "; var at xi = 0.5 and 0.7:",
    format(P$var_fn(list(mu = c(0, 0), sigma = c(1, 1), xi = c(0.5, 0.7)),
                    list())), "\n")
# brms's posterior_epred_gen_extreme_value(), the formula as brms writes it
be <- function(mu, sigma, xi) mu + sigma * (gamma(1 - xi) - 1) / xi
cat("brms's epred formula at xi = 0, 1e-8, 1, 1.5 (mu 0, sigma 1):",
    format(suppressWarnings(be(0, 1, c(0, 1e-8, 1, 1.5))), digits = 6),
    "\n")
cat("frmtmb's mean_fn at the same points:                     ",
    format(P$mean_fn(list(mu = rep(0, 4), sigma = rep(1, 4),
                          xi = c(0, 1e-8, 1, 1.5)), list()), digits = 6),
    "\n")
cat("brms's own function on a prep (brms:::posterior_epred_gen_extreme_value):",
    format(brms:::posterior_epred_gen_extreme_value(list(dpars = list(
      mu = rep(0, 4), sigma = rep(1, 4), xi = c(0, 1e-8, 1, 1.5)))),
      digits = 6), "\n")
lb <- grid$mu - grid$sigma
ub <- grid$mu + 2 * grid$sigma
cdf <- function(q, m, s, k) {
  x <- (q - m) / s
  if (k == 0) exp(-exp(-x)) else exp(-pmax(1 + k * x, 0)^(-1 / k))
}
tm_ref <- mapply(function(m, s, k, a, b) {
  q_int(function(p) qf(p, m, s, k), cdf(a, m, s, k), cdf(b, m, s, k)) /
    (cdf(b, m, s, k) - cdf(a, m, s, k))
}, grid$mu, grid$sigma, grid$xi, lb, ub)
cat("trunc_mean_fn on (mu - sigma, mu + 2 sigma): max rel err",
    format(max(rel(P$trunc_mean_fn(as.list(grid), list(), lb, ub), tm_ref)),
           digits = 3), "\n")

cat("\n=== discrete_weibull\n")
fam <- discrete_weibull()
P <- post(fam)
grid <- expand.grid(mu = c(0.05, 0.5, 0.9, 0.99, 0.999),
                    k = c(0.3, 0.6, 1, 2, 5))
# reference: the direct sum to 1e7 terms, plus the tail beyond it as the
# integral from 1e7 + 1/2 (the midpoint rule, not the Euler-Maclaurin form
# frmtmb uses), which at these parameters is below 1e-13 of the sum except
# where the column `tail_share` says otherwise
ref <- mapply(function(mu, k) {
  y <- seq_len(1e7)
  S <- exp(log(mu) * y^k)
  L <- -log(mu)
  tl <- exp(lgamma(1 / k) + stats::pgamma(L * (1e7 + 0.5)^k, 1 / k,
                                          lower.tail = FALSE, log.p = TRUE) -
              log(k) - log(L) / k)
  c(mean = sum(S) + tl, m2 = sum((2 * y - 1) * S), tail_share =
      tl / (sum(S) + tl))
}, grid$mu, grid$k)
mf <- P$mean_fn(list(mu = grid$mu, shape = grid$k), list())
vf <- P$var_fn(list(mu = grid$mu, shape = grid$k), list())
ok <- ref["tail_share", ] < 1e-13
cat("mean_fn vs the 1e7-term sum: max rel err",
    format(max(rel(mf[ok], ref["mean", ok])), digits = 3), "over", sum(ok),
    "of", nrow(grid), "points (the rest have a tail beyond 1e7 terms)\n")
okv <- ok & ref["tail_share", ] < 1e-16
v_ref <- ref["m2", ] - ref["mean", ]^2
cat("var_fn vs the 1e7-term sum: max rel err",
    format(max(rel(vf[okv], v_ref[okv])), digits = 3), "over", sum(okv),
    "points\n")
cat("where the tail beyond 1e7 terms matters (mu, shape, frmtmb mean,",
    "sum + tail integral):\n")
for (i in which(!ok)) {
  cat(sprintf("  mu = %g shape = %g: %.12g vs %.12g, rel err %.2e\n",
              grid$mu[i], grid$k[i], mf[i], ref["mean", i],
              rel(mf[i], ref["mean", i])))
}
# brms's mean: the sum stopped at the first term below 0.001
bm <- vapply(seq_len(nrow(grid)), function(i) {
  suppressWarnings(brms:::mean_discrete_weibull(grid$mu[i], grid$k[i]))
}, 0)
r <- rel(bm, mf)
cat("brms:::mean_discrete_weibull() against frmtmb's mean: rel diff",
    "median", format(median(r), digits = 3), "max", format(max(r), digits = 3),
    "at mu =", grid$mu[which.max(r)], "shape =", grid$k[which.max(r)],
    "(brms", format(bm[which.max(r)], digits = 8), "frmtmb",
    format(mf[which.max(r)], digits = 8), ")\n")
cat("brms warns that its approximation failed at",
    sum(vapply(seq_len(nrow(grid)), function(i) {
      w <- FALSE
      withCallingHandlers(brms:::mean_discrete_weibull(grid$mu[i],
                                                       grid$k[i]),
                          warning = function(e) {
                            w <<- TRUE
                            invokeRestart("muffleWarning")
                          })
      w
    }, NA)), "of", nrow(grid), "points\n")
lb <- c(1, 2, 3)
ub <- c(4, Inf, 10)
g3 <- list(mu = c(0.6, 0.9, 0.95), shape = c(0.8, 1.3, 0.5))
tm <- P$trunc_mean_fn(g3, list(), lb, ub)
tm_ref <- vapply(1:3, function(i) {
  y <- lb[i]:min(ub[i], 1e6)
  p <- exp(log(g3$mu[i]) * y^g3$shape[i]) -
    exp(log(g3$mu[i]) * (y + 1)^g3$shape[i])
  sum(y * p) / sum(p)
}, 0)
cat("trunc_mean_fn (lb, ub) = (1, 4), (2, Inf), (3, 10) vs direct sums:",
    "max rel err", format(max(rel(tm, tm_ref)), digits = 3), "\n")

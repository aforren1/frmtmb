# Probe 2: does the three-alternative image series tape, is its gradient
# right, and what does it cost?
#
# Run from the worktree root:
#   Rscript dev/gddm-nchoice/probe2-tape.R
#
# Compare every number against probe0's shipped 1-D gddm: 8.56 s to build
# the tape for one condition at dt = 0.01, ny = 201, t_max = 2, and 6.7 ms
# per gradient over 2000 trials.

source("dev/gddm-nchoice/common.R")
source("dev/gddm-nchoice/images.R")
suppressPackageStartupMessages(library(numDeriv))

tl_all <- gd_tiles(400L)
trim <- function(m) list(A = tl_all$A[seq_len(m)],
                         w = tl_all$w[seq_len(m), , drop = FALSE],
                         sgn = tl_all$sgn[seq_len(m)], geom = tl_all$geom)

# --- where the truncated series stops being usable ------------------------
#
# The image sum is an alternating sum of Gaussians. At a decision time t the
# images that still matter reach out to about 6 sqrt(t), so the number of
# tiles a horizon needs grows like the AREA of that disc over the area of
# the triangle: linearly in t / r^2. Past the truncation the sum is not
# merely inaccurate, it is garbage, because the terms that would have
# cancelled are missing.
cat("-- usable horizon against tile count, zero drift, start at centre --\n")
cat("reported: the largest t at which the m-tile sum is within 1e-8 of the\n")
cat("400-tile sum, in units of the dimensionless time t / r^2\n\n")
cat(sprintf("%8s", "tiles"))
rs <- c(0.6, 0.866, 1.6)
for (r in rs) cat(sprintf(" %14s", sprintf("r = %.3f", r)))
cat(sprintf(" %10s\n", "t / r^2"))
tv <- exp(seq(log(0.02), log(30), length.out = 240))
for (m in c(22L, 46L, 76L, 118L, 166L, 250L, 340L)) {
  tm <- trim(m)
  lim <- numeric(length(rs))
  for (q in seq_along(rs)) {
    r <- rs[q]
    a <- vapply(tv, function(t)
      gd_img_density(t, c(0, 0), c(0, 0), r, tm, 1L), numeric(1))
    b <- vapply(tv, function(t)
      gd_img_density(t, c(0, 0), c(0, 0), r, tl_all, 1L), numeric(1))
    ok <- abs(a - b) <= 1e-8 * pmax(abs(b), 1e-12)
    lim[q] <- if (all(ok)) Inf else tv[which(!ok)[1L]]
  }
  cat(sprintf("%8d", m))
  for (q in seq_along(rs)) cat(sprintf(" %14.3f", lim[q]))
  cat(sprintf(" %10.2f\n", mean(lim / rs^2)))
}

# --- the likelihood, taped ------------------------------------------------
#
# Free parameters: two drift contrasts (the third is minus their sum, the
# model only sees relative evidence), the threshold c, a two-vector start
# bias, and the non-decision time. Six, against wiener()'s four.
mk_data <- function(ntrial, seed = 11L) {
  set.seed(seed)
  mu <- c(0.9, 0.1, -1.0); cc <- 0.7071; ndt <- 0.25
  tlq <- tl_all
  a <- as.numeric(t(gd_basis(3L)) %*% mu); r <- sqrt(1.5) * cc
  tg <- seq(1e-4, 12, length.out = 30000L)
  f <- lapply(1:3, function(k) gd_img_density(tg, a, c(0, 0), r, tlq, k))
  P <- vapply(f, function(v) sum(diff(tg) * (v[-1L] + v[-length(v)]) / 2),
              numeric(1))
  ch <- sample.int(3L, ntrial, TRUE, prob = P / sum(P))
  rt <- numeric(ntrial)
  for (k in 1:3) {
    idx <- which(ch == k)
    cdf <- c(0, cumsum(diff(tg) * (f[[k]][-1L] + f[[k]][-30000L]) / 2)) / P[k]
    rt[idx] <- stats::approx(cdf, tg, xout = stats::runif(length(idx)))$y
  }
  list(choice = ch, rt = rt + ndt,
       truth = list(mu = mu, cc = cc, ndt = ndt, a = a, r = r))
}

mk_obj <- function(d, tlm) {
  who <- lapply(1:3, function(k) which(d$choice == k))
  U <- gd_basis(3L)
  f <- function(p) {
    a <- as.numeric(t(U) %*% c(p$mu12, -sum(p$mu12)))
    r <- sqrt(1.5) * exp(p$logc)
    nd <- p$ndt
    ll <- 0
    for (k in 1:3) {
      if (!length(who[[k]])) next
      ll <- ll + sum(log(gd_img_density(d$rt[who[[k]]] - nd, a, p$z0, r,
                                        tlm, k) + 1e-300))
    }
    -ll
  }
  RTMB::MakeADFun(f, list(mu12 = c(0.9, 0.1), logc = log(0.7071),
                          z0 = c(0, 0), ndt = 0.25), silent = TRUE)
}

cat("\n-- gradient against numDeriv, 200 trials, 118 tiles --\n")
d200 <- mk_data(200L)
o <- mk_obj(d200, trim(118L))
x <- o$par + c(0.05, -0.03, 0.02, 0.01, -0.01, 0.005)
g1 <- o$gr(x)
g2 <- numDeriv::grad(function(v) o$fn(v), x,
                     method = "Richardson")
cat(sprintf("%8s %16s %16s %12s\n", "par", "RTMB", "numDeriv", "rel"))
for (i in seq_along(x))
  cat(sprintf("%8s %16.8f %16.8f %12.2e\n", names(o$par)[i], g1[i], g2[i],
              abs(g1[i] - g2[i]) / max(abs(g2[i]), 1e-8)))

cat("\n-- Hessian: does the second order tape too? --\n")
h1 <- o$he(x)
h2 <- numDeriv::jacobian(function(v) o$gr(v), x)
cat(sprintf("max |RTMB he - numDeriv jac(gr)| / max|he|: %.2e\n",
            max(abs(h1 - h2)) / max(abs(h1))))

cat("\n-- cost: tape build and one gradient --\n")
cat(sprintf("%8s %8s %12s %12s %12s %14s\n",
            "trials", "tiles", "build s", "fn s", "grad s", "nll"))
for (nt in c(500L, 2000L, 8000L)) {
  dd <- mk_data(nt)
  for (m in c(46L, 118L, 250L)) {
    tb <- gd_time(mk_obj(dd, trim(m)))
    ob <- tb$value
    tf <- gd_time(ob$fn(ob$par), reps = 5L)
    tg2 <- gd_time(ob$gr(ob$par), reps = 5L)
    cat(sprintf("%8d %8d %12.3f %12.5f %12.5f %14.3f\n",
                nt, m, tb$sec, tf$sec, tg2$sec, tf$value))
  }
}

cat("\n-- does it recover the truth? 4000 trials, 118 tiles --\n")
d4 <- mk_data(4000L, seed = 21L)
o4 <- mk_obj(d4, trim(118L))
t0 <- proc.time()[["elapsed"]]
fit <- stats::nlminb(o4$par, o4$fn, o4$gr,
                     control = list(iter.max = 400L, eval.max = 600L))
el <- proc.time()[["elapsed"]] - t0
sdr <- try(RTMB::sdreport(o4), silent = TRUE)
se <- if (inherits(sdr, "try-error")) rep(NA_real_, 6) else sqrt(diag(sdr$cov.fixed))
tr <- d4$truth
truth <- c(tr$mu[1], tr$mu[2], log(tr$cc), 0, 0, tr$ndt)
cat(sprintf("nlminb: %.1f s, %d iterations, conv %d, objective %.4f\n",
            el, fit$iterations, fit$convergence, fit$objective))
cat(sprintf("%8s %12s %12s %10s %8s\n", "par", "estimate", "truth", "se", "z"))
for (i in seq_along(fit$par))
  cat(sprintf("%8s %12.5f %12.5f %10.5f %8.2f\n", names(o4$par)[i],
              fit$par[i], truth[i], se[i], (fit$par[i] - truth[i]) / se[i]))

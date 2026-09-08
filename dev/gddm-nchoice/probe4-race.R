# Probe 4: does route 1, the race families already in the package, already
# cover what a genuine multi-alternative diffusion would buy?
#
# Run from the worktree root:
#   Rscript dev/gddm-nchoice/probe4-race.R
#
# Three measurements.
#
# 1. A QUALITATIVE separation. The triangle model is a RELATIVE-evidence
#    model: only mu_k - mean(mu) enters, so adding the same constant to
#    every drift changes nothing at all. A race with absolute thresholds
#    speeds up. A design that raises overall evidence while holding
#    relative evidence fixed therefore tells the two apart outright, and
#    neither family can imitate the other there.
#
# 2. Fit rdm(3) and lba(3) to a large sample from the triangle model and
#    read the maximized log-likelihood against the triangle model's own.
#    For a correctly specified model the gap to the truth is about p/2, so
#    anything much larger is misspecification, not sampling noise.
#
# 3. Say what that misfit costs in observable terms: the choice proportions
#    and the response time quantiles the best-fitting race implies, against
#    the truth, in units of the Monte Carlo error of a realistic experiment.

source("dev/gddm-nchoice/common.R")
source("dev/gddm-nchoice/images.R")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})

tl <- gd_tiles(250L)
tl_fit <- local({ m <- 118L
  list(A = tl$A[seq_len(m)], w = tl$w[seq_len(m), , drop = FALSE],
       sgn = tl$sgn[seq_len(m)], geom = tl$geom) })
U <- gd_basis(3L)
tg <- seq(1e-4, 12, length.out = 24000L)
trapc <- function(v) c(0, cumsum(diff(tg) * (v[-1L] + v[-length(v)]) / 2))

tri_pieces <- function(mu, cc, z0 = c(0, 0)) {
  a <- as.numeric(t(U) %*% mu); r <- sqrt(1.5) * cc
  f <- lapply(1:3, function(k) gd_img_density(tg, a, z0, r, tl, k))
  P <- vapply(f, function(v) trapc(v)[length(tg)], numeric(1))
  list(f = f, P = P, a = a, r = r)
}

tri_sim <- function(n, mu, cc, ndt, z0 = c(0, 0), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  s <- tri_pieces(mu, cc, z0)
  ch <- sample.int(3L, n, TRUE, prob = s$P / sum(s$P))
  rt <- numeric(n)
  for (k in 1:3) {
    idx <- which(ch == k)
    if (!length(idx)) next
    cdf <- trapc(s$f[[k]]) / s$P[k]
    rt[idx] <- stats::approx(cdf, tg, xout = stats::runif(length(idx)))$y
  }
  data.frame(rt = rt + ndt, choice = ch)
}

qtl <- function(d, p = c(0.1, 0.5, 0.9))
  do.call(rbind, lapply(1:3, function(k) {
    v <- d$rt[d$choice == k]
    c(P = length(v) / nrow(d), stats::quantile(v, p, names = FALSE))
  }))

# --- 1. the qualitative separation ---------------------------------------
cat("-- 1. a common shift in every drift: relative-evidence model vs race --\n")
mu <- c(0.9, 0.1, -1.0); cc <- 0.7071; ndt <- 0.25
for (sh in c(0, 1, 3)) {
  d <- tri_sim(60000L, mu + sh, cc, ndt, seed = 5L)
  q <- qtl(d)
  cat(sprintf("triangle, drift + %d: P = %.4f %.4f %.4f   median rt = %.4f\n",
              sh, q[1, 1], q[2, 1], q[3, 1], stats::median(d$rt)))
}
for (sh in c(0, 1, 3)) {
  v <- c(0.9, 0.6, 0.3) + sh
  d <- rdm_simulate(60000L, v = v, A = 0.5, k = 0.7, ndt = ndt)
  q <- qtl(d)
  cat(sprintf("rdm,      drift + %d: P = %.4f %.4f %.4f   median rt = %.4f\n",
              sh, q[1, 1], q[2, 1], q[3, 1], stats::median(d$rt)))
}

# --- 2. fit the races to triangle data ------------------------------------
cat("\n-- 2. maximized log-likelihood, 20000 trials from the triangle --\n")
N <- 20000L
dat <- tri_sim(N, mu, cc, ndt, seed = 9L)
cat(sprintf("data: P = %s, median rt %.4f\n",
            paste(sprintf("%.4f", table(dat$choice) / N), collapse = " "),
            stats::median(dat$rt)))

tri_obj <- function(d) {
  who <- lapply(1:3, function(k) which(d$choice == k))
  f <- function(p) {
    a <- as.numeric(t(U) %*% c(p$mu12, -sum(p$mu12)))
    r <- sqrt(1.5) * exp(p$logc)
    ll <- 0
    for (k in 1:3) ll <- ll + sum(log(
      gd_img_density(d$rt[who[[k]]] - p$ndt, a, p$z0, r, tl_fit, k) + 1e-300))
    -ll
  }
  RTMB::MakeADFun(f, list(mu12 = mu[1:2], logc = log(cc), z0 = c(0, 0),
                          ndt = ndt), silent = TRUE)
}
o <- tri_obj(dat)
ll_truth <- -o$fn(o$par)
t0 <- proc.time()[["elapsed"]]
ft <- stats::nlminb(o$par, o$fn, o$gr,
                    control = list(iter.max = 300L, eval.max = 500L))
cat(sprintf("triangle at truth  logLik %12.3f\n", ll_truth))
cat(sprintf("triangle at MLE    logLik %12.3f  (%d par, %.0f s, conv %d)\n",
            -ft$objective, length(ft$par), proc.time()[["elapsed"]] - t0,
            ft$convergence))

for (fam in c("rdm", "lba")) {
  t0 <- proc.time()[["elapsed"]]
  fit <- try(frm(bf(rt | vint(choice) ~ 1),
                 family = if (fam == "rdm") rdm(3) else lba(3), data = dat),
             silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat(sprintf("%-8s FAILED: %s\n", fam, conditionMessage(attr(fit, "condition"))))
    next
  }
  lv <- as.numeric(stats::logLik(fit))
  cat(sprintf("%-8s at MLE    logLik %12.3f  (%d par, %.0f s)  gap to triangle %.1f\n",
              fam, lv, length(stats::coef(fit)),
              proc.time()[["elapsed"]] - t0, -ft$objective - lv))
  assign(paste0("fit_", fam), fit)
}

# --- 3. what the misfit looks like ---------------------------------------
cat("\n-- 3. the best-fitting race against the truth, in observable terms --\n")
cat("z is the discrepancy in units of the Monte Carlo error of a 500-trial\n")
cat("experiment, so |z| below about 2 would go unnoticed in one subject\n\n")
qd <- qtl(dat)
for (fam in c("rdm", "lba")) {
  fo <- get0(paste0("fit_", fam))
  if (is.null(fo)) next
  sm <- stats::simulate(fo, nsim = 1L, seed = 3L)
  ch <- sm[[1L]]
  # simulate() returns response times; the choice is held fixed by the
  # family's own draw, so read the implied distribution by resampling the
  # fitted family directly instead
  cf <- stats::coef(fo)
  cat(sprintf("%s fitted coefficients: %s\n", fam,
              paste(sprintf("%s=%.3f", names(cf), cf), collapse = " ")))
}
lnk <- function(fo) {
  cf <- stats::coef(fo)
  v <- exp(cf[grep("^v", names(cf))])
  list(v = as.numeric(v), A = exp(cf[["A_Intercept"]]),
       k = exp(cf[["k_Intercept"]]),
       ndt = frmtmb::single_response(fo)$family$links$ndt$linkinv(cf[["ndt_Intercept"]]))
}
if (!is.null(get0("fit_rdm"))) {
  pr <- try(lnk(fit_rdm), silent = TRUE)
  if (!inherits(pr, "try-error")) {
    sim <- rdm_simulate(200000L, v = pr$v, A = pr$A, k = pr$k, ndt = pr$ndt)
    qs <- qtl(sim)
    cat(sprintf("\n%8s %10s %10s %10s %10s\n", "", "truth", "rdm fit",
                "se(500)", "z"))
    lab <- c("P", "q10", "q50", "q90")
    for (k in 1:3) for (j in 1:4) {
      tv <- qd[k, j]; fv <- qs[k, j]
      se <- if (j == 1) sqrt(tv * (1 - tv) / 500) else
        stats::sd(dat$rt[dat$choice == k]) / sqrt(500 * qd[k, 1]) *
        c(1.71, 1.25, 1.71)[j - 1]
      cat(sprintf("k=%d %-4s %10.4f %10.4f %10.4f %10.2f\n",
                  k, lab[j], tv, fv, se, (fv - tv) / se))
    }
  }
}

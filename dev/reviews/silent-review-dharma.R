# Review of item 650. Usage: Rscript silent-review-dharma.R <part> [base]
# part: hand | power | trunc | edge
args <- commandArgs(TRUE)
part <- args[1]
arm <- if (length(args) > 1) args[2] else "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
ks_p <- function(r) suppressWarnings(stats::ks.test(r, "punif")$p.value)

# right censoring at a row-varying time drawn to give censored share ~ rate
sim_right <- function(seed, n = 300, rate = 0.5, skew = FALSE) {
  set.seed(seed)
  x <- rnorm(n)
  e <- if (skew) (rexp(n) - 1) else rnorm(n)
  ys <- 1 + 0.5 * x + e
  # censoring time c = q + U(-1, 1), q the (1 - rate) quantile of ys
  q <- stats::quantile(ys, 1 - rate)
  ct <- q + runif(n, -1, 1)
  cen <- as.numeric(ys > ct)
  data.frame(x = x, y = ifelse(cen == 1, ct, ys), cen = cen)
}

if (part == "hand") {
  # Dunn & Smyth by hand from the ML estimates vs dharma_residuals()
  for (s in 1:5) {
    d <- sim_right(s, rate = 0.4)
    fit <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = d)
    cf <- coef(fit)
    mu <- fitted(fit)
    if (is.matrix(mu) || is.data.frame(mu)) mu <- mu[, 1]
    sig <- as.numeric(sigma(fit))
    F <- pnorm(d$y, mu, sig)
    set.seed(1000 + s)
    uh <- ifelse(d$cen == 1, runif(nrow(d), F, 1), F)
    ud <- dharma_residuals(fit, nsim = 2000, seed = s)$scaledResiduals
    # per-row: expected value of the dharma residual on censored rows is
    # (Fhat(y) + 1)/2; compare its rank interval with the exact one
    sims <- as.matrix(simulate(fit, nsim = 2000, seed = s))
    lo_emp <- rowMeans(sims <= d$y)
    i <- d$cen == 1
    cat(sprintf(paste0("HAND seed %d cens %.2f sigma %.4f | KS unif hand p %.3f",
                       " dharma p %.3f | 2-sample hand vs dharma p %.3f |",
                       " max|Femp - F| censored rows %.4f",
                       " | corr(ud, uh) uncens %.4f\n"),
                s, mean(d$cen), sig, ks_p(uh), ks_p(ud),
                suppressWarnings(stats::ks.test(uh, ud)$p.value),
                max(abs(lo_emp[i] - F[i])), cor(ud[!i], uh[!i])))
  }
}

if (part == "power") {
  for (rate in c(0.1, 0.5, 0.9)) {
    pc <- pw <- numeric(20)
    for (s in 1:20) {
      d <- sim_right(s, rate = rate)
      f <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = d)
      pc[s] <- ks_p(dharma_residuals(f, nsim = 250, seed = s)$scaledResiduals)
      dw <- sim_right(s, rate = rate, skew = TRUE)
      fw <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = dw)
      pw[s] <- ks_p(dharma_residuals(fw, nsim = 250, seed = s)$scaledResiduals)
    }
    cat(sprintf(paste0("POWER %s rate %.1f: correct rejected p<0.01 %d/20",
                       " p<0.05 %d/20 (median p %.3f); wrong rejected",
                       " p<0.01 %d/20 p<0.05 %d/20 (median p %.3g)\n"),
                arm, rate, sum(pc < 0.01), sum(pc < 0.05), median(pc),
                sum(pw < 0.01), sum(pw < 0.05), median(pw)))
  }
}

if (part == "trunc") {
  # cens() with trunc(): gaussian truncated below at 0, right-censored
  pr <- pl <- pp <- numeric(20)
  for (s in 1:20) {
    set.seed(s)
    n <- 600
    x <- rnorm(n)
    ys <- 1 + 0.5 * x + rnorm(n)
    keep <- ys > 0
    x <- x[keep][1:300]; ys <- ys[keep][1:300]
    ct <- runif(300, 0.8, 2.5)
    cen <- as.numeric(ys > ct)
    d <- data.frame(x = x, y = ifelse(cen == 1, ct, ys), cen = cen)
    f <- frm(bf(y | cens(cen) + trunc(lb = 0) ~ x), family = gaussian(),
             data = d)
    pr[s] <- ks_p(dharma_residuals(f, nsim = 250, seed = s)$scaledResiduals)
    # left censoring under an upper truncation
    set.seed(100 + s)
    x <- rnorm(n)
    ys <- 1 + 0.5 * x + rnorm(n)
    keep <- ys < 2
    x <- x[keep][1:300]; ys <- ys[keep][1:300]
    lc <- runif(300, -0.5, 1)
    cen <- -as.numeric(ys < lc)
    d <- data.frame(x = x, y = ifelse(cen == -1, lc, ys), cen = cen)
    f <- frm(bf(y | cens(cen) + trunc(ub = 2) ~ x), family = gaussian(),
             data = d)
    pl[s] <- ks_p(dharma_residuals(f, nsim = 250, seed = s)$scaledResiduals)
    # poisson, interval (y, y2] and left censoring on the lattice
    set.seed(200 + s)
    x <- rnorm(300)
    ys <- rpois(300, exp(1.2 + 0.4 * x))
    cen <- numeric(300); y <- ys; y2 <- ys
    il <- ys <= 1 & runif(300) < 0.7
    cen[il] <- -1; y[il] <- 1
    # cells (2, 4], (4, 6], ...: coarsening at random, ys >= 3 only
    ii <- ys >= 3 & runif(300) < 0.4
    # interval cells (2j, 2j + 2]: y = 2 * floor((ys - 1) / 2) ... event
    # y < Y <= y2 with y2 = y + 2 on a partition of the lattice
    lo <- 2 * floor((ys - 1) / 2) + 1 - 1
    cen[ii] <- 2; y[ii] <- lo[ii]; y2[ii] <- lo[ii] + 2
    d <- data.frame(x = x, y = y, y2 = y2, cen = cen)
    f <- frm(bf(y | cens(cen, y2) ~ x), family = poisson(), data = d)
    pp[s] <- ks_p(dharma_residuals(f, nsim = 250, seed = s)$scaledResiduals)
  }
  cat(sprintf("TRUNC %s right+lb: p<0.01 %d/20 median %.3f\n", arm,
              sum(pr < 0.01), median(pr)))
  cat(sprintf("TRUNC %s left+ub: p<0.01 %d/20 median %.3f\n", arm,
              sum(pl < 0.01), median(pl)))
  cat(sprintf("TRUNC %s poisson left+interval: p<0.01 %d/20 median %.3f\n",
              arm, sum(pp < 0.01), median(pp)))
}

if (part == "edge") {
  d <- sim_right(3, rate = 0.3)
  d$w <- runif(nrow(d), 0.5, 2)
  r <- try(dharma_residuals(frm(bf(y | cens(cen) + weights(w) ~ x),
                                family = gaussian(), data = d),
                            nsim = 100, seed = 1))
  cat("EDGE weights:", if (inherits(r, "try-error")) "ERROR" else "ran", "\n")
  # mi() response with cens()
  d2 <- d
  d2$y[c(5, 7)] <- NA
  d2$cen[5] <- 1
  f2 <- try(frm(bf(y | mi() + cens(cen) ~ x), family = gaussian(), data = d2))
  if (!inherits(f2, "try-error")) {
    r2 <- try(dharma_residuals(f2, nsim = 100, seed = 1))
    cat("EDGE mi+cens:", if (inherits(r2, "try-error")) "ERROR" else
      paste("ran, n =", length(r2$scaledResiduals)), "\n")
  } else cat("EDGE mi+cens: fit refused\n")
  # NA cens with na.action exclude
  d3 <- d
  d3$x[4] <- NA
  f3 <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = d3)
  r3 <- try(dharma_residuals(f3, nsim = 100, seed = 1))
  cat("EDGE NA covariate row:", if (inherits(r3, "try-error")) "ERROR" else
    paste("ran, n =", length(r3$scaledResiduals)), "\n")
  # rotation through ...
  r4 <- try(dharma_residuals(frm(bf(y | cens(cen) ~ x), family = gaussian(),
                                 data = d), nsim = 100, seed = 1,
                             rotation = "estimated"))
  cat("EDGE rotation:", if (inherits(r4, "try-error")) "ERROR" else "ran", "\n")
  # gaussian left censoring by a fixed point (type I) and right at fixed point
  set.seed(5)
  x <- rnorm(300); ys <- 1 + 0.5 * x + rnorm(300)
  cen <- as.numeric(ys > 1.5)
  dd <- data.frame(x = x, y = pmin(ys, 1.5), cen = cen)
  ff <- frm(bf(y | cens(cen) ~ x), family = gaussian(), data = dd)
  rr <- dharma_residuals(ff, nsim = 250, seed = 2)
  cat(sprintf("EDGE type-I right at 1.5: cens %.2f KS p %.3f\n", mean(cen),
              ks_p(rr$scaledResiduals)))
  # a lognormal (positive support), right-censored
  pl <- numeric(10)
  for (s in 1:10) {
    set.seed(s)
    x <- rnorm(300); ys <- exp(0.3 + 0.4 * x + 0.6 * rnorm(300))
    ct <- runif(300, 0.5, 3)
    cen <- as.numeric(ys > ct)
    dl <- data.frame(x = x, y = ifelse(cen == 1, ct, ys), cen = cen)
    fl <- frm(bf(y | cens(cen) ~ x), family = lognormal(), data = dl)
    pl[s] <- ks_p(dharma_residuals(fl, nsim = 250, seed = s)$scaledResiduals)
  }
  cat(sprintf("EDGE lognormal right: p<0.01 %d/10 median %.3f\n",
              sum(pl < 0.01), median(pl)))
}

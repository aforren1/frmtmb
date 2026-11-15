# Lane nlstop2: reruns of two unclear items on the arm's build.
#   3359: an active upper bound on the top coefficient of a degree-8 raw
#         polynomial (dev/illse2-p4bound.R's design, seeds 1 to 3).
#   2818: n = 2e4, y ~ s(x, k = 100) + (1 | g), gaussian
#         (dev/seflags-rev-cost.R smooth_big, set.seed(1)).
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-rerun.R <3359|2818>
source("dev/nlstop2-common.R")
what <- commandArgs(TRUE)[1]
if (what == "3359") {
  for (seed in 1:3) {
    deg <- 8
    set.seed(seed)
    d <- data.frame(x = runif(200, 1, 2))
    d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
    xn <- paste0("x", 1:deg)
    for (k in 1:deg) d[[xn[k]]] <- d$x^k
    top <- xn[deg]
    fo <- as.formula(paste("y ~", paste0(xn, collapse = " + ")))
    l0 <- lm(fo, data = d)
    ub <- coef(l0)[[top]] - sqrt(diag(vcov(l0)))[[top]]
    t0 <- proc.time()[["elapsed"]]
    r <- tt(frm(fo, data = d,
                prior = set_prior("", class = "b", coef = top, ub = ub)))
    t1 <- proc.time()[["elapsed"]] - t0
    d$off <- ub * d[[top]]
    fo1 <- as.formula(paste("y ~", paste0(xn[-deg], collapse = " + "),
                            "+ offset(off)"))
    l1 <- lm(fo1, data = d)
    cat(sprintf("3359 deg 8 seed %d: logLik - constrained lm %+.4g code %s t %.2f warn %s\n",
                seed, llv(r) - as.numeric(logLik(l1)),
                if (iserr(r)) "ERR" else r$opt$convergence, t1, wsum(r)))
    if (iserr(r)) cat("  ", r$msg, "\n")
  }
}
if (what == "2818") {
  set.seed(1)
  n <- 2e4
  d <- data.frame(x = runif(n), g = factor(sample(100, n, TRUE)))
  d$y <- sin(8 * d$x) + rnorm(100, 0, 0.3)[d$g] + rnorm(n, 0, 0.5)
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(y ~ s(x, k = 100) + (1 | g), family = gaussian(), data = d))
  t1 <- proc.time()[["elapsed"]] - t0
  cat(sprintf("2818: ll %.6f code %s t %.1f warn %s\n", llv(r),
              if (iserr(r)) "ERR" else r$opt$convergence, t1, wsum(r)))
  if (iserr(r)) cat("  ", r$msg, "\n")
  for (w in attr(r, "w")) cat("  W:", substr(w, 1, 200), "\n")
  if (!iserr(r) && requireNamespace("mgcv", quietly = TRUE)) {
    m <- mgcv::gamm(y ~ s(x, k = 100), random = list(g = ~1), data = d,
                    method = "ML")
    cat(sprintf("2818: mgcv::gamm ML logLik %.6f\n", as.numeric(logLik(m$lme))))
  }
}
cat("DONE\n")

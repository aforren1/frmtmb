# Lane nlstop2: the new assertions of test-fit-restart.R on 20 draws of
# the data a few ulps apart (x * (1 + 2e-15 k), k = 1 to 20).
#   pole:  binomial mm seeds 11 and 44 (fresh rows): the default fit's
#          logLik against the start at K = -6, and the message class
#   user:  bernoulli mm 55 from Vm -2, K -2.5: the pole warning
#   explin: a * (exp(b * x) - 1), n 1000, seed 2, start (2, 0.5): the
#          limit warning
#   huber: huber() y ~ x + (1 | g), seed 4, scale 1: the default fit
#          against the best of 4 jittered starts
# Usage: NLSTOP2_ARM=before|after Rscript dev/nlstop2-ulp.R <what>
source("dev/nlstop2-common.R")
fresh_defs()
what <- commandArgs(TRUE)[1]
ks <- 1:20
if (what == "pole") {
  for (seed in c(11, 44)) {
    d <- gen("binomial", "mm", seed * 1000 + 1)
    fm <- form("binomial", "mm")
    ref <- frm(fm, family = F$binomial, data = d,
               start = list(beta = c(-0.1, -6)))
    for (k in ks) {
      dk <- d
      dk$xp <- dk$xp * (1 + 2e-15 * k)
      r <- tt(frm(fm, family = F$binomial, data = dk))
      cat(sprintf("pole seed %d k %d: ref - ll %.3g K %.5g code %s w %s m %s\n",
                  seed, k, llv(ref) - llv(r), r$estimates$beta[2],
                  r$opt$convergence, wsum(r), msum(r)))
    }
  }
}
if (what == "user") {
  d <- gen("bernoulli", "mm", 55 * 1000 + 1)
  fm <- form("bernoulli", "mm")
  for (k in ks) {
    dk <- d
    dk$xp <- dk$xp * (1 + 2e-15 * k)
    r <- tt(frm(fm, family = F$bernoulli, data = dk,
                start = list(beta = c(-2, -2.5))))
    cat(sprintf("user k %d: ll %.6f K %.5g code %s w %s\n", k, llv(r),
                r$estimates$beta[2], r$opt$convergence, wsum(r)))
  }
}
if (what == "explin") {
  set.seed(2)
  d <- data.frame(x = runif(1000, 0, 3))
  d$y <- 2 * d$x + rnorm(1000, 0, 0.5)
  for (k in ks) {
    dk <- d
    dk$x <- dk$x * (1 + 2e-15 * k)
    r <- tt(frm(bf(y ~ a * (exp(b * x) - 1), a + b ~ 1, nl = TRUE),
                data = dk, start = list(beta = c(2, 0.5))))
    cat(sprintf("explin k %d: ll %.6f code %s w %s\n", k, llv(r),
                r$opt$convergence, wsum(r)))
  }
}
if (what == "huber") {
  set.seed(4)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  b <- rnorm(20)[d$g]
  d$y <- 1 + 0.5 * d$x + 0.5 * b + rt(100, 3)
  for (k in ks) {
    dk <- d
    dk$y <- dk$y * (1 + 2e-15 * k)
    r <- tt(frm(y ~ x + (1 | g), family = huber(), data = dk))
    best <- llv(r)
    set.seed(4 * 31 + 7)
    for (i in 1:4) {
      st <- lapply(r$estimates[c("beta", "betad", "theta")], function(v) {
        unname(v + rnorm(length(v), 0, 0.3) * pmax(1, abs(v)))
      })
      best <- max(best, llv(tt(frm(y ~ x + (1 | g), family = huber(),
                                   data = dk, start = st))), na.rm = TRUE)
    }
    cat(sprintf("huber k %d: quad %s best - ll %.3g w %s\n", k,
                isTRUE(r$quadrature), best - llv(r), wsum(r)))
  }
}
cat("DONE\n")

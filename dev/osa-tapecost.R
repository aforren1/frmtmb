# Lane osa, fix round: what one tape evaluation costs, by the CDF pair it
# carries. 500 rows in 25 groups, linear predictor u[g] + 0.5, counts
# from the matching rng (seed 1). Each arm is the density plus both log
# tails at every row, as the one-step tape carries them. Time of 300
# gradient-and-Hessian evaluations at moving u, minimum of 3 rounds.
#   Rscript dev/osa-tapecost.R
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(frmtmb) })
set.seed(1)
n <- 500; G <- 25; k <- 10; shp <- 2.5
g <- rep(seq_len(G), each = n / G)
e <- rnorm(G, 0, 0.7)[g]
yp <- rpois(n, exp(0.5 + e))
yb <- rbinom(n, k, plogis(0.5 + e))
y1 <- rbinom(n, 1, plogis(0.5 + e))
yn <- rnbinom(n, size = shp, mu = exp(0.5 + e))
u0 <- rnorm(G, 0, 0.1)
pq <- function(u) {
  eta <- u[g] + 0.5
  list(p = exp(-logspace_add(0 * eta, -eta)),
       q = exp(-logspace_add(0 * eta, eta)))
}
arms <- list(
  pois_density = function(u) sum(dpois(yp, exp(u[g] + 0.5), log = TRUE)),
  pois_summed = function(u) {
    mu <- exp(u[g] + 0.5)
    lp <- frmtmb:::poisson_cdf_log(yp, mu)
    sum(dpois(yp, mu, log = TRUE)) + sum(lp$lF) + sum(lp$lS)
  },
  pois_pgamma = function(u) {
    mu <- exp(u[g] + 0.5)
    S <- pgamma(mu, shape = yp + 1)
    sum(dpois(yp, mu, log = TRUE)) + sum(log(S)) + sum(log(1 - S))
  },
  binom_density = function(u) {
    sum(dbinom(yb, k, plogis(u[g] + 0.5), log = TRUE))
  },
  binom_ibeta = function(u) {
    z <- pq(u)
    lp <- frmtmb:::binomial_cdf_log(yb, k, z)
    sum(dbinom(yb, k, z$p, log = TRUE)) + sum(lp$lF) + sum(lp$lS)
  },
  binom_pbeta = function(u) {
    z <- pq(u)
    qc <- ifelse(yb >= k, 0, yb)
    sum(dbinom(yb, k, z$p, log = TRUE)) +
      sum(log(pbeta(z$q, k - qc, qc + 1))) +
      sum(log(pbeta(z$p, qc + 1, k - qc)))
  },
  bern_ibeta = function(u) {
    z <- pq(u)
    lp <- frmtmb:::binomial_cdf_log(y1, 1, z)
    sum(dbinom(y1, 1, z$p, log = TRUE)) + sum(lp$lF) + sum(lp$lS)
  },
  bern_closed = function(u) {
    z <- pq(u)
    sum(dbinom(y1, 1, z$p, log = TRUE)) + sum(log(z$q) * (y1 == 0)) +
      sum(log(z$p) * (y1 == 0))
  },
  nb_density = function(u) {
    mu <- exp(u[g] + 0.5)
    sum(dnbinom2(yn, mu, mu + mu^2 / shp, log = TRUE))
  },
  nb_ibeta = function(u) {
    lmu <- u[g] + 0.5
    lp <- frmtmb:::nbinom_cdf_log(yn, lmu, log(shp) + 0 * lmu)
    mu <- exp(lmu)
    sum(dnbinom2(yn, mu, mu + mu^2 / shp, log = TRUE)) + sum(lp$lF) +
      sum(lp$lS)
  },
  nb_pbeta = function(u) {
    lmu <- u[g] + 0.5
    pr <- exp(-logspace_add(0 * lmu, lmu - log(shp)))
    p1 <- exp(-logspace_add(0 * lmu, log(shp) - lmu))
    mu <- exp(lmu)
    sum(dnbinom2(yn, mu, mu + mu^2 / shp, log = TRUE)) +
      sum(log(pbeta(pr, shp, yn + 1))) + sum(log(pbeta(p1, yn + 1, shp)))
  })
cat("BEGIN osa-tapecost\n")
for (nm in names(arms)) {
  tp <- MakeTape(arms[[nm]], u0)
  h <- tp$jacfun()
  hh <- h$jacfun()
  best <- Inf
  for (r in 1:3) {
    t0 <- proc.time()[["elapsed"]]
    for (j in 1:300) {
      uj <- u0 + j * 1e-6
      h(uj)
      hh(uj)
    }
    best <- min(best, proc.time()[["elapsed"]] - t0)
  }
  cat(sprintf("%-14s 300 gradients and Hessians %7.3fs\n", nm, best))
}
cat("END osa-tapecost\n")

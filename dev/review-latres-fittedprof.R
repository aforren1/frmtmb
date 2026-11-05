# Attack 6: fitted() took 11 s at 2000 rows against a 0.3 s fit. Which step, and
# is it the lane's? Control: the same inner dimension as a random effect of one
# level per row, (1 | id), and as ar1(ft + 0 | g) with 20 groups of 100.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
for (n in c(500, 1000, 2000, 4000)) {
  set.seed(1)
  d <- data.frame(t = seq_len(n), g = factor(1), id = factor(seq_len(n)), x = rnorm(n))
  e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
  d$y <- rpois(n, exp(.5 + .3 * d$x + e))
  fl <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  fc <- suppressWarnings(frm(y ~ x + (1 | id), family = poisson(), data = d))
  t1 <- system.time(a <- fitted(fl))[[3]]
  t2 <- system.time(b <- fitted(fc))[[3]]
  t3 <- system.time(frm_linpred(fl))[[3]]
  t4 <- system.time(frm_linpred(fl, se.fit = TRUE))[[3]]
  t5 <- system.time(frm_linpred(fc, se.fit = TRUE))[[3]]
  cat(sprintf("n %5d  fitted: latent %7.2f  (1|id) %7.2f   linpred no se: latent %6.2f   linpred se: latent %7.2f (1|id) %7.2f\n", n, t1, t2, t3, t4, t5))
}
n <- 2000
set.seed(1)
d <- data.frame(t = seq_len(n), g = factor(1), id = factor(seq_len(n)), x = rnorm(n))
e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
d$y <- rpois(n, exp(.5 + .3 * d$x + e))
fl <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
tmp <- tempfile()
Rprof(tmp, interval = 0.05, memory.profiling = TRUE)
a <- fitted(fl)
Rprof(NULL)
s <- summaryRprof(tmp, memory = "both")
print(head(s$by.self, 8))
print(head(s$by.total[order(-s$by.total$total.time), c("total.time", "self.time", "mem.total")], 25))

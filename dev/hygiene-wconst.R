# Item 131: margins of test-weights-constant.R's two assertions over 20
# perturbations of the response by a few ulps (y * (1 + 2e-15 * k)),
# k = 0..19, seed 31 as in the test.
#   Rscript dev/hygiene-wconst.R
.libPaths(c(if (nzchar(Sys.getenv("HYG_LIB"))) Sys.getenv("HYG_LIB"),
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(31)
n <- 80L
d0 <- data.frame(x = stats::rnorm(n), z = stats::runif(n))
d0$y <- 1 + 0.6 * d0$x - 0.4 * d0$z + stats::rnorm(n, 0, 0.7)
out <- NULL
for (k in 0:19) {
  d <- d0
  d$y <- d$y * (1 + 2e-15 * k)
  d$w2 <- 2
  d$w3 <- 3
  f1 <- frm(bf(y ~ x + z), family = gaussian(), data = d)
  f2 <- frm(bf(y | weights(w2) ~ x + z), family = gaussian(), data = d)
  f1s <- frm(bf(y ~ x + z, sigma ~ x), family = gaussian(), data = d)
  f3 <- frm(bf(y | weights(w3) ~ x + z, sigma ~ x), family = gaussian(),
            data = d)
  r1 <- max(abs(f2$opt$par - f1$opt$par) / sqrt(diag(vcov(f1, full = TRUE))))
  r2 <- max(abs(f3$opt$par - f1s$opt$par) /
              sqrt(diag(vcov(f1s, full = TRUE))))
  l1 <- as.numeric(logLik(f1))
  l1s <- as.numeric(logLik(f1s))
  q1 <- abs(as.numeric(logLik(f2)) - 2 * l1) / abs(l1)
  q2 <- abs(as.numeric(logLik(f3)) - 3 * l1s) / abs(l1s)
  out <- rbind(out, data.frame(k, r1, r2, q1, q2))
}
print(out, digits = 3)
cat(sprintf("MAX est/se %.3g (bound 1e-3)  MAX rel logLik %.3g (bound 1e-8)\n",
            max(out$r1, out$r2), max(out$q1, out$q2)))
cat(sprintf("logLik at k=0: unweighted %.4f, w=2 %.4f\n", l1,
            as.numeric(logLik(f2))))

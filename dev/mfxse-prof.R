# Where one get_predict() call's time goes on a poisson (1 + x | g) fit at 500
# levels and 20000 rows (data seed 7411): Rprof over 300 calls with a moved b.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
set.seed(7411)
J <- 500L; m <- 40L; n <- J * m
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(J), each = m)))
d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(J, 0, 0.8)[d$g]))
fit <- quiet(frm(y ~ x + (1 + x | g), family = poisson(), data = d))
b0 <- fit$estimates$b
pf <- file.path(mfx_root, "dev/mfxse-prof.out")
Rprof(pf, interval = 0.002)
t0 <- proc.time()[["elapsed"]]
for (i in 1:300) {
  f <- fit
  f$estimates$b[i] <- b0[i] + 1e-6
  get_predict(f, newdata = d)
}
el <- proc.time()[["elapsed"]] - t0
Rprof(NULL)
cat(sprintf("%.4f s per get_predict\n", el / 300))
s <- summaryRprof(pf)
print(head(s$by.total, 30))
print(head(s$by.self, 15))

# The dense joint covariance behind fitted(): time and peak memory by n, latent
# against a random effect of one level per row (the control), one fit per process.
# Usage: Rscript review-latres-dense.R <n> <latent|re>
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
a <- commandArgs(TRUE)
n <- as.integer(a[[1]])
set.seed(1)
d <- data.frame(t = seq_len(n), g = factor(1), id = factor(seq_len(n)), x = rnorm(n))
d$y <- rpois(n, exp(.5 + .3 * d$x + as.numeric(arima.sim(list(ar = .5), n, sd = .5))))
peak <- function() as.numeric(system(sprintf('powershell -NoProfile -Command "(Get-Process -Id %d).PeakWorkingSet64"', Sys.getpid()), intern = TRUE)[1]) / 1e9
f <- suppressWarnings(if (a[[2]] == "latent") frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d) else frm(y ~ x + (1 | id), family = poisson(), data = d))
p0 <- peak()
t1 <- system.time(v <- fitted(f))[[3]]
p1 <- peak()
tmp <- tempfile()
saveRDS(f, tmp)
cat(sprintf("n %d %s fitted %.1f s  peak before %.2f GB after %.2f GB  rds %.0f MB\n", n, a[[2]], t1, p0, p1, file.size(tmp) / 1e6))

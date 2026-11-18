# Lane nlband: which ending does the ridge test of test-fit-restart.R
# ("a body on a ridge to its linear limit says so", bernoulli mm seed 55)
# take on this R? The test has one more assertion when the fit warns
# frmtmb_nl_limit. Run with the reference R and an OpenBLAS copy.
# Rscript dev/nlband-obridge.R
source("dev/nlband-common.R")
cat("R.home:", R.home(), "\n")
set.seed(55 * 1000 + 1)
n <- 300
d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                xq = round(runif(n, 0, 5), 3),
                xb = round(runif(n, 0, 3000), 1),
                g = factor(sample(1:15, n, TRUE)))
d$ug <- rnorm(15, 0, 0.4)[d$g]
lat <- 3 * d$xp / (1 + d$xp)
d$y <- rbinom(n, 1, plogis(lat - 2))
r <- tt(frm(bf(y ~ Vm * xp / (K + xp), Vm + K ~ 1, nl = TRUE),
            family = bernoulli(), data = d))
cat("RIDGE ll", format(llv(r), digits = 12), "code", r$opt$convergence,
    "warn", wsum(r), "msg", msum(r), "\n")

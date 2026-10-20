# Lane adhess: is the point sdreport() reads (last.par.best) the fit's
# own estimate on the x20 exp fit, and how far does the gap move the
# exact covariance? Seed 1 of dev/adhess-nls.R.
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
set.seed(1)
d <- data.frame(x = 20 + runif(200))
d$y <- 2 * exp(0.8 * (d$x - 20)) + rnorm(200, 0, 0.3)
f <- frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE), data = d,
         start = list(beta = c(2 * exp(-16), 0.8)))
p <- f$opt$par
q <- ns$sdr_outer_point(f$obj)
cat("opt$par - last.par.best:", format(p - q, digits = 3), "\n")
H1 <- f$obj$he(p); H2 <- f$obj$he(q)
print(cov_err(solve(H2), solve(H1)))
print(cov_err(suppressWarnings(vcov(f, full = TRUE)), solve(H2)))

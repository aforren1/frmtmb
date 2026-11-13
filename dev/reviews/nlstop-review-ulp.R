# Lane nlstop, item 4850: the linear fit under normal(0, 1000) on the slope
# at y about 1e6 (seed 6201), the response moved by a few ulps over 20
# draws: does every draw reach the mode near the prior's location, with
# the message?
# Review copy of dev/nlstop-4850-ulp.R with the seed as argument. Usage: Rscript dev/reviews/nlstop-review-ulp.R <seed>
source("dev/nlstop-common.R")
set.seed(as.integer(commandArgs(TRUE)[1]))
n <- 200
d0 <- data.frame(x = rnorm(n))
d0$y <- 0.5 * d0$x + rnorm(n)
d0$yl <- 1e6 * (1 + 0.05 * d0$x) + rnorm(n, 0, 1e3)
pr <- c(set_prior("normal(0, 1000)", class = "b"),
        set_prior("normal(1e6, 1e6)", class = "Intercept"))
v <- numeric(0); nm <- 0L
for (k in 1:20) {
  d <- d0
  d$yl <- d$yl * (1 + 2e-15 * k)
  r <- tt(frm(yl ~ x, data = d, prior = pr))
  v <- c(v, llp(r))
  nm <- nm + any(attr(r, "mc") == "frmtmb_prior_multistart")
  cat(sprintf("draw %d: llp %.6f slope %.2f code %s warn %s\n", k, llp(r),
              r$estimates$beta[2], r$opt$convergence, wsum(r)))
}
cat(sprintf("20 draws: llp from %.6f to %.6f; message on %d\n", min(v),
            max(v), nm))

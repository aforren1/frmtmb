# Review of lane nlstop: the case prior_conflict_escape() excludes, a
# prior on the centered intercept in conflict with the data. y ~ x at y
# about 1e6 (seed 6201, the 4850 data), normal(0, psd) on the Intercept.
# The default fit against a start at the prior's location (Intercept 0,
# the slope at its data value, log sigma at log(sd(y)) + 6, i.e. inflated);
# then rstan::optimizing() on brms 2.23.0's program from both points.
# Usage: Rscript dev/reviews/nlstop-review-intercept.R
source("dev/nlstop-common.R")
suppressPackageStartupMessages({library(brms); library(rstan)})
source("tests/testthat/helper-brms.R")
set.seed(6201); n <- 200
d <- data.frame(x = rnorm(n))
d$y <- 0.5 * d$x + rnorm(n)
d$yl <- 1e6 * (1 + 0.05 * d$x) + rnorm(n, 0, 1e3)
cat("ROWHEAD psd arm llp code warn msg b_Int b sigma\n")
for (psd in c(1e3, 3e3, 1e4, 2e4, 3e4, 1e5)) {
  pr <- set_prior(sprintf("normal(0, %g)", psd), class = "Intercept")
  r0 <- tt(frm(yl ~ x, data = d, prior = pr))
  # start: beta = (centered or not, frmtmb's own Intercept), slope;
  # betad = log sigma
  st <- list(beta = c(0, 5e4), betad = log(1e6))
  r1 <- tt(frm(yl ~ x, data = d, prior = pr, start = st))
  for (a in list(list("default", r0), list("start0", r1))) {
    r <- a[[2]]
    if (iserr(r)) { cat("ROW", psd, a[[1]], "ERR", r$msg, "\n"); next }
    fe <- fixef(r)[, "Estimate"]
    cat(sprintf("ROW %g %s %.6f %s %s %s %.6g %.6g %.6g\n", psd, a[[1]],
                llp(r), r$opt$convergence, wsum(r),
                if (length(attr(r, "mc"))) paste(unique(attr(r, "mc")), collapse = "+") else "-",
                fe[1], fe[2], exp(r$estimates$betad[1])))
  }
  if (psd %in% c(1e4, 3e4)) {
    bpr <- c(pr, set_prior("", class = "sigma"))
    mod <- brms_stan_model(make_stancode(yl ~ x, data = d, prior = bpr))
    sdat <- make_standata(yl ~ x, data = d, prior = bpr)
    for (a in list(list("default", r0), list("start0", r1))) {
      r <- a[[2]]
      fe <- fixef(r)[, "Estimate"]
      init <- list(Intercept = fe[[1]] + fe[[2]] * mean(d$x),
                   b = array(fe[[2]]),
                   sigma = exp(r$estimates$betad[[1]]))
      o <- suppressWarnings(optimizing(mod, data = sdat, init = init,
        algorithm = "LBFGS", iter = 20000, as_vector = FALSE,
        verbose = FALSE))
      cat(sprintf("  stan from %s: lp %.6f rc %d b_Int %.6g b %.6g sigma %.6g\n",
                  a[[1]], o$value, o$return_code,
                  o$par$Intercept - o$par$b * mean(d$x), o$par$b[1],
                  o$par$sigma))
    }
  }
  flush.console()
}
cat("DONE\n")

# Lane gpfield fix round (copy of dev/reviews/gpfield-review-cost.R): cost of predict() on an exact gp()
# with n positions, lane against base, in separate processes.
# Usage: Rscript gpfield-review-cost.R <lane|base> <fam> <n> <ndraws>
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
a <- commandArgs(TRUE)
arm <- a[1]; fam <- a[2]; n <- as.integer(a[3]); nd_ <- as.integer(a[4])
gf_libs(arm)
set.seed(9300 + n)
d <- data.frame(x = sort(round(stats::runif(n, 0, 50), 4)))
eta <- sin(d$x / 3)
d$y <- switch(fam, gaussian = 2 * eta + stats::rnorm(n, 0, 0.3),
              poisson = stats::rpois(n, exp(1 + eta)))
ff <- switch(fam, gaussian = gaussian(), poisson = poisson())
t0 <- proc.time()[["elapsed"]]
q <- gf_fit(bf(y ~ gp(x)), data = d, family = ff)
fit <- q$fit
cat(sprintf("COST %s %s n %d fit wall %.1f conv %d ll %.6f warn [%s]\n", arm,
            fam, n, proc.time()[["elapsed"]] - t0, fit$opt$convergence,
            as.numeric(logLik(fit)), paste(unique(sub(":.*", "", q$warn)),
                                           collapse = ";")))
newd <- data.frame(x = c(d$x[c(10, 200)], 25.05, 51, 55, 80))
for (k in 1:2) {
  set.seed(1)
  t0 <- proc.time()[["elapsed"]]
  p <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd_)
  w <- proc.time()[["elapsed"]] - t0
  cat(sprintf("COST %s %s n %d predict ndraws %d wall %.2f per draw %.4f | sd %s\n",
              arm, fam, n, nd_, w, w / nd_,
              paste(round(apply(p, 2, sd), 3), collapse = " ")))
}
set.seed(1)
t0 <- proc.time()[["elapsed"]]
p <- predict(fit, summary = FALSE, ndraws = nd_)
w <- proc.time()[["elapsed"]] - t0
cat(sprintf("COST %s %s n %d predict(no newdata) ndraws %d wall %.2f\n", arm,
            fam, n, nd_, w))
cat("mem", gc()[2, 6], "MB max\n")

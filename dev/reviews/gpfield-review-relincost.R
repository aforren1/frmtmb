# Review of lane gpfield: the cost of gp_hsgp_relin() and gp_mean_relin()
# on larger fits, lane against base (separate processes).
# Usage: Rscript gpfield-review-relincost.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
arm <- commandArgs(TRUE)[1]
gf_libs(arm)
set.seed(9900)
n <- 5000
d <- data.frame(x = runif(n, 0, 10), x2 = runif(n, 0, 10))
d$y <- sin(d$x) + 0.5 * cos(d$x2 / 2) + rnorm(n, 0, 0.5)
d$yp <- rpois(n, exp(0.5 + 0.5 * sin(d$x)))
rows <- list(
  h1g = list(bf(y ~ gp(x, k = 30)), gaussian()),
  h2g = list(bf(y ~ gp(x, x2, k = 12)), gaussian()),
  h1p = list(bf(yp ~ gp(x, k = 30)), poisson()))
for (nm in names(rows)) {
  t0 <- proc.time()[["elapsed"]]
  q <- gf_fit(rows[[nm]][[1]], data = d, family = rows[[nm]][[2]])
  w <- proc.time()[["elapsed"]] - t0
  f <- q$fit
  cat(sprintf("RCOST %s %s n %d wall %.2f ll %.6f conv %d gp2 %s warn [%s]\n",
              arm, nm, n, w, as.numeric(logLik(f)), f$opt$convergence,
              paste(format(f$opt$gp_second %||% NA, digits = 3),
                    collapse = ","),
              paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
}
# the mean's reading on an exact poisson field, n 400
set.seed(9901)
n <- 400
d <- data.frame(x = sort(runif(n, 0, 20)))
d$y <- rpois(n, exp(1 + sin(d$x)))
t0 <- proc.time()[["elapsed"]]
q <- gf_fit(bf(y ~ gp(x)), data = d, family = poisson())
f <- q$fit
cat(sprintf("RCOST %s exact-pois n %d wall %.2f ll %.6f conv %d gp2 %s warn [%s]\n",
            arm, n, proc.time()[["elapsed"]] - t0, as.numeric(logLik(f)),
            f$opt$convergence,
            paste(format(f$opt$gp_second %||% NA, digits = 3), collapse = ","),
            paste(unique(sub(":.*", "", q$warn)), collapse = ";")))

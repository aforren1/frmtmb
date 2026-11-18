# Lane alapse, fix round 1: what the quantile search costs. y ~ x,
# gaussian noise, quantile estimated, n = 100, 1000, 5000 (seed 7);
# elapsed time of frm() (minimum of 3) on the base build and the lane,
# and the lane's linear programs (ald_cert$search$evals).
#   Rscript dev/alapse-searchcost.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
for (n in c(100, 1000, 5000)) {
  set.seed(7)
  d <- data.frame(x = rnorm(n))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  tt <- numeric(3)
  for (k in 1:3) {
    tt[k] <- system.time(f <- suppressWarnings(frm(y ~ x,
      family = asym_laplace(), data = d)))[["elapsed"]]
  }
  s <- f$opt$ald_cert$search
  cat(sprintf("n %5d | frm() %.2f s (min of 3) | code %d | LPs %s | global %s\n",
              n, min(tt), f$opt$convergence, format(s$evals %||% NA),
              format(isTRUE(f$opt$ald_cert$global))))
}

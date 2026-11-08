# Lane arstep: is the student step row's objective path dependent? The
# Laplace objective at ONE outer point, from different inner starts: zeros,
# the default fit's mode, and the mode of the best reference fit (start sd
# 1.5, lscale 0.1, dev/arstep-diag3.R). Data as dev/review-ncgp-r3-ng.R
# (seed 74000 + 31 seed + n). Args: arm case n seed
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]; cs <- args[2]; n <- as.integer(args[3]); seed <- as.integer(args[4])
rev_libs(arm)
f <- function(x) 2 * (x > 5) - 1
set.seed(74000 + seed * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
d <- data.frame(x = x, y = f(x) + 0.02 * rt(n, 5))
r0 <- rev_fit(y ~ gp(x), d, family = student())
r1 <- rev_fit(y ~ gp(x), d, family = student(), start = list(theta = c(log(1.5), log(0.1))))
cat(sprintf("IN default ll %.4f; reference ll %.4f\n", as.numeric(logLik(r0$fit)), as.numeric(logLik(r1$fit))))
for (who in c("default", "reference")) {
  fit <- if (who == "default") r0$fit else r1$fit
  obj <- fit$obj
  env <- obj$env
  rnd <- env$random
  p <- fit$opt$par
  modes <- list(zeros = numeric(length(rnd)), default = r0$fit$obj$env$last.par.best[rnd],
                reference = r1$fit$obj$env$last.par.best[rnd])
  for (nm in names(modes)) {
    lp <- env$last.par.best
    lp[rnd] <- modes[[nm]]
    env$last.par.best <- lp
    env$last.par <- lp
    v <- obj$fn(p)
    cat(sprintf("IN at the %s optimum, inner start %-9s -> objective %.4f\n", who, nm, v))
  }
}
cat("DONE\n")

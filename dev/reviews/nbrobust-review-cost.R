# Review of lane nbrobust, check 5: time per objective and per gradient
# of a negbinomial GLMM, y ~ x + (1 | g), on one build. Different data
# from the lane's: seed 21, n = 3000, 60 groups, shape 1.5, mean about 6.
#   Rscript dev/reviews/nbrobust-review-cost.R <base|lib> <rep>
args <- commandArgs(TRUE)
arm <- args[1]
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/wt-nbrobust-lib", base)
          else base)
suppressPackageStartupMessages(library(frmtmb))
set.seed(21)
n <- 3000
g <- sample(60, n, TRUE)
x <- stats::rnorm(n)
u <- stats::rnorm(60, sd = 0.5)
d <- data.frame(x = x, g = factor(g),
                y = stats::rnbinom(n, size = 1.5, mu = exp(1.5 + 0.3 * x + u[g])))
t0 <- proc.time()[["elapsed"]]
f <- frm(y ~ x + (1 | g), family = negbinomial(), data = d)
tfit <- proc.time()[["elapsed"]] - t0
p <- f$opt$par
timeit <- function(fun, minsec = 2) {
  best <- Inf
  for (r in 1:5) {
    k <- 0L
    t0 <- proc.time()[["elapsed"]]
    repeat {
      fun()
      k <- k + 1L
      el <- proc.time()[["elapsed"]] - t0
      if (el >= minsec) break
    }
    best <- min(best, el / k)
  }
  best
}
ctl <- timeit(function() sum(sort(stats::runif(2e4))))
tfn <- timeit(function() f$obj$fn(p))
tgr <- timeit(function() f$obj$gr(p))
cat(sprintf(paste("RESULT %s rep %s logLik %.8f lsh %.6f fit %.2f s",
                  "control %.4e fn %.4e gr %.4e fn/ctl %.4f gr/ctl %.4f\n"),
            arm, args[2], as.numeric(logLik(f)), p[["betad"]], tfit, ctl,
            tfn, tgr, tfn / ctl, tgr / ctl))

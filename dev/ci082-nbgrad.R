# Lane ci082: the gradient of RTMB::dnbinom_robust() in log(var - mu)
# at a large shape, against the exact derivative.
#   Rscript dev/ci082-nbgrad.R
# With size n = exp(-s) * mu^2 ... the log-density's derivative in the
# shape is tiny there (the poisson limit), so any large value is noise.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
cat("RTMB", as.character(packageVersion("RTMB")), " TMB",
    as.character(packageVersion("TMB")), "\n")
set.seed(2)
x <- stats::rnorm(200)
y <- stats::rpois(200, exp(1 + 0.4 * x))
lmu <- 1.0544184774 + 0.4005673777 * x
f <- function(p) {
  sum(RTMB::dnbinom_robust(y, lmu, 2 * lmu - p[1], log = TRUE))
}
tp <- RTMB::MakeTape(f, 20)
gtp <- tp$jacfun()
# exact: d/d(log shape) of sum log NB(y | mu, shape), with
# shape = exp(ls): shape * sum(digamma(y + shape) - digamma(shape) +
# log(shape / (shape + mu)) + 1 - (y + shape) / (shape + mu)), the
# digamma difference summed as sum_{j < y} 1 / (shape + j) and the log
# term as -log1p(mu / shape), so nothing cancels
exact <- function(ls) {
  s <- exp(ls)
  mu <- exp(lmu)
  dg <- vapply(y, function(yi) if (yi == 0) 0 else
    sum(1 / (s + seq_len(yi) - 1)), 0)
  s * sum(dg - log1p(mu / s) + (mu - y) / (s + mu))
}
for (ls in c(10, 14, 18, 20, 22, 24, 24.4791351853, 26, 28, 30, 35)) {
  cat(sprintf("log shape %8.4f  tape grad % .4e  exact % .4e\n", ls,
              gtp(ls), exact(ls)))
}

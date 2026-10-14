# Lane illse: a nonlinear formula carrying a degree-8 raw polynomial
# (each coefficient its own nonlinear parameter), where the design's
# coordinates do not apply: the SE check must keep its verdict.
#   ILLSE_ARM=base|lane Rscript dev/illse-nl.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
for (seed in 1:3) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  f <- bf(y ~ b0 + b1 * x + b2 * x^2 + b3 * x^3 + b4 * x^4 + b5 * x^5 +
            b6 * x^6 + b7 * x^7 + b8 * x^8,
          b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 ~ 1, nl = TRUE)
  r <- cap(frm(f, data = d))
  m <- r$v
  cat(sprintf("seed %d code %d lost %s root %s\n  %s\n", seed,
              m$opt$convergence, lost_of(m),
              !is.null(frmtmb:::sdr_of(m)$cov_fixed_root),
              paste(substr(r$w, 1, 160), collapse = "\n  ")))
}

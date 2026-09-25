# Punch round 2, item 3: Rmpfr interval masses off the review's set,
# which starts at u1 = 0.05 and |v| a = 24. These go to u1 = 1e-3 and
# |v| a = 80, with the review's 1200-bit eigenfunction tail
# (dev/phase3b-review2/r2-mpfr-lib.R).
# Usage: Rscript dev/phase3b-interval-ref-small.R
# Output: dev/phase3b-log/interval-small.rds
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(Rmpfr))
source("dev/phase3b-review2/r2-mpfr-lib.R")
set.seed(2027)
n <- 240
pts <- data.frame(v = runif(n, -20, 20), a = exp(runif(n, log(0.3), log(4))),
                  w = runif(n, 0.05, 0.95),
                  u1 = exp(runif(n, log(1e-3), log(0.05))),
                  r = 1 + exp(runif(n, log(1e-3), log(3))))
pts$t1 <- pts$u1 * pts$a^2
pts$t2 <- pts$t1 * pts$r
pts$ref <- vapply(seq_len(n), function(i) {
  as.numeric(log(tail_eigen(pts$t1[i], pts$v[i], pts$a[i], pts$w[i]) -
                 tail_eigen(pts$t2[i], pts$v[i], pts$a[i], pts$w[i])))
}, 0)
pts$refF1 <- vapply(seq_len(n), function(i) {
  as.numeric(log(p_lower(pts$v[i], pts$a[i], pts$w[i]) -
                 tail_eigen(pts$t1[i], pts$v[i], pts$a[i], pts$w[i])))
}, 0)
pts$mass_over_F1 <- exp(pts$ref - pts$refF1)
saveRDS(pts, "dev/phase3b-log/interval-small.rds")
cat("points:", n, " finite refs:", sum(is.finite(pts$ref)), "\n")

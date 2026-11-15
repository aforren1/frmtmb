# Short check, item 4: cost of one fit per process (peak working set is
# that fit's). The lane's design (dev/arspat2-cost2.R, from the review):
# rook lattice r x c, row-standardized W (sar error), M =
# solve(crossprod(I - 0.4 W)) (fcor), (1 | h) 8 levels, student() with
# gamma(2, 0.1) on nu, data seed 11, missing seed m.
#   Rscript dev/reviews/arspat2-short-cost.R <rows> <cols> <m> <fcor|sar>
# Base when SHORT_LIB="".
a <- commandArgs(TRUE)
nr <- as.integer(a[1]); nc <- as.integer(a[2]); m <- as.integer(a[3])
ty <- a[4]
lib <- Sys.getenv("SHORT_LIB", "C:/Users/adf44/source/r/rellib-r22")
.libPaths(c(if (nzchar(lib)) lib else "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
arm <- if (nzchar(lib)) "lane" else "base"
g <- expand.grid(r = seq_len(nr), c = seq_len(nc))
n <- nrow(g)
Aj <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- Aj / rowSums(Aj)
set.seed(11)
d <- data.frame(x = rnorm(n), h = factor(rep(1:8, length.out = n)))
A0 <- diag(n) - 0.4 * W
d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(8)[d$h] +
                         rnorm(n) * sqrt(4 / rchisq(1, 4))))
M <- if (ty == "fcor") solve(crossprod(A0))
set.seed(m)
d$y[sort(sample(n, m))] <- NA
pr <- set_prior("gamma(2, 0.1)", class = "nu")
fo <- if (ty == "fcor") y | mi() ~ x + (1 | h) + fcor(M) else
  y | mi() ~ x + (1 | h) + sar(W, type = "error")
d2 <- if (ty == "fcor") list(M = M) else list(W = W)
peak <- function() ps::ps_memory_info()[["peak_wset"]] / 2^30
gc()
p0 <- peak()
t0 <- proc.time()
fit <- tryCatch(suppressWarnings(suppressMessages(
  frm(fo, data = d, data2 = d2, family = student(), prior = pr))),
  error = function(e) e)
tt <- proc.time() - t0
if (inherits(fit, "error")) {
  cat(sprintf("COST %s n %d m %d %s: ERROR %s, %.2f s, peak %.2f GB\n", arm,
              n, m, ty, conditionMessage(fit), tt[["elapsed"]], peak()))
} else {
  gg <- fit$obj$gr(fit$opt$par)
  cat(sprintf("COST %s n %d m %d %s: fit %.2f s, logLik %.4f, conv %d, grad finite %s, peak %.2f GB (before fit %.2f)\n",
              arm, n, m, ty, tt[["elapsed"]], as.numeric(logLik(fit)),
              fit$opt$convergence, all(is.finite(gg)), peak(), p0))
}

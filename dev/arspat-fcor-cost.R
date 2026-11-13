# Lane arspat, backlog 4120: what fcor()'s pointwise log-density costs per
# draw, on the new route (diag(M^-1) read off the block) and the old one
# (K = Linv / sigma and K * K formed per draw), in ONE process on the lane
# build. The old route is the same build with `diag_inv` taken off the
# block, which sends autocor_seam_loo_parts() down the generic jacobian
# path it took before. Peak memory from gc()'s "max used" after a reset;
# time from blocks of repeated calls grown past 1.2 s, minimum of 5
# rounds, the two arms interleaved; the control is the new route timed
# against itself. Output: dev/arspat-out/fcor-cost.txt
.libPaths(c("C:/Users/adf44/source/r/wt-arspat-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
n <- as.integer(Sys.getenv("ARSPAT_N", "1500"))
d <- data.frame(x = rnorm(n))
d$y <- 1 + d$x + rnorm(n)
M <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
fit <- frm(y ~ x + fcor(M), data = d, data2 = list(M = M))
old <- fit
old$frame$autocor$y$aux$diag_inv <- NULL
dpv <- frmtmb::eval_dpars(fit)
l_new <- frmtmb::autocor_row_loglik(fit, dpv, "y")
l_old <- frmtmb::autocor_row_loglik(old, dpv, "y")
cat(sprintf("n = %d; max relative difference of the two routes %.2e\n", n,
            max(abs(l_new - l_old)) / max(abs(l_old))))
peak <- function(f) {
  invisible(gc(reset = TRUE))
  base <- gc()[2, 6]
  f()
  gc()[2, 6] - base
}
cat(sprintf("peak Vcells memory over the call, Mb: new %.1f old %.1f\n",
            peak(function() frmtmb::autocor_row_loglik(fit, dpv, "y")),
            peak(function() frmtmb::autocor_row_loglik(old, dpv, "y"))))
reps <- 1L
repeat {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) frmtmb::autocor_row_loglik(fit, dpv, "y")
  if (proc.time()[["elapsed"]] - t0 > 1.2) break
  reps <- reps * 2L
}
arm <- function(f) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) frmtmb::autocor_row_loglik(f, dpv, "y")
  (proc.time()[["elapsed"]] - t0) / reps
}
tt <- replicate(5, c(new = arm(fit), old = arm(old), control = arm(fit)))
m <- apply(tt, 1, min)
cat(sprintf("seconds per draw (min of 5 rounds of %d): new %.4f old %.4f control %.4f\n",
            reps, m[["new"]], m[["old"]], m[["control"]]))
cat(sprintf("ratio old / new %.2f; control / new %.2f\n",
            m[["old"]] / m[["new"]], m[["control"]] / m[["new"]]))

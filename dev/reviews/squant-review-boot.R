# Review of lane squant, item 4: the parametric bootstrap MEDIUM on 50
# data sets (homo design of dev/squant-se.R: n 400, y = 1 + 2 x + N(0, 1),
# quantile 0.9, smoothing 0.25), 100 draws each, seeds 20261200 + i.
# Records the slope, its Hessian, sandwich and bootstrap SEs, and the
# coverage of each 95% Wald interval for the true slope 2.
# Run (worktree root):
#   Rscript dev/reviews/squant-review-boot.R > dev/reviews/squant-review-boot.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
t0 <- proc.time()
rows <- list()
for (i in 1:50) {
  set.seed(20261200 + i)
  n <- 400
  dd <- data.frame(x = runif(n), row = seq_len(n))
  dd$y <- 1 + 2 * dd$x + rnorm(n)
  fit <- frm(bf(y ~ x, quantile = 0.9), family = quantile_elf(), data = dd)
  b <- suppressWarnings(frm_bootstrap(fit, nsim = 100, seed = i))
  bs <- stats::sd(b$t[b$converged, 2])
  rows[[i]] <- data.frame(
    i = i, conv = fit$opt$convergence, b1 = fixef(fit)[2, "Estimate"],
    hess = sqrt(vcov(fit)[2, 2]),
    sand = sqrt(vcov_cluster(fit, ~ row, "CR0")[2, 2]), boot = bs,
    nconv = sum(b$converged))
  cat(sprintf("data set %d: slope %.4f | SE hessian %.4f sandwich %.4f boot %.4f (%d of 100 refits converged)\n",
              i, rows[[i]]$b1, rows[[i]]$hess, rows[[i]]$sand, bs,
              rows[[i]]$nconv))
}
r <- do.call(rbind, rows)
mcsd <- stats::sd(r$b1)
cover <- function(se) mean(abs(r$b1 - 2) <= 1.959964 * se)
cat(sprintf(paste0("SUMMARY 50 data sets: MC sd of the slope %.4f (lane's 300",
                   " sets: 0.2871) | mean SE / MC sd: hessian %.3f sandwich",
                   " %.3f bootstrap %.3f | mean SE / 0.2871: bootstrap %.3f |",
                   " coverage hessian %.2f sandwich %.2f bootstrap %.2f |",
                   " bootstrap / hessian median %.3f | refits converged %d of",
                   " %d | %.0f s\n"),
            mcsd, mean(r$hess) / mcsd, mean(r$sand) / mcsd,
            mean(r$boot) / mcsd, mean(r$boot) / 0.2871, cover(r$hess),
            cover(r$sand), cover(r$boot), stats::median(r$boot / r$hess),
            sum(r$nconv), 50 * 100, (proc.time() - t0)[3]))
cat("DONE\n")

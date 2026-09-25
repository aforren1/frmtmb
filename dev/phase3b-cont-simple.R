# Item 3.5: is lambda's Wald interval right on the simplest design?
# One subject, n trials, mu 0.8, bs 1.4, ndt 0.25, bias 0.5, 5 percent
# contaminants uniform on [lo, hi] with a coin-flip boundary.
# Usage: Rscript dev/phase3b-cont-simple.R <n> <nrep> <lo> <hi> <maxndt|none>
.libPaths(c(Sys.getenv("P3B_LIB", "C:/Users/adf44/source/r/phase3b-lib2"),
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
a <- commandArgs(trailingOnly = TRUE)
n <- as.integer(a[[1]]); nrep <- as.integer(a[[2]])
lo <- as.numeric(a[[3]]); hi <- as.numeric(a[[4]])
mx <- if (a[[5]] == "none") NULL else as.numeric(a[[5]])
res <- NULL
for (seed in seq_len(nrep)) {
  set.seed(seed)
  d <- ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.25, bias = 0.5)
  hit <- runif(n) < 0.05
  d$rt[hit] <- runif(sum(hit), lo, hi)
  d$upper[hit] <- rbinom(sum(hit), 1, 0.5)
  fit <- suppressWarnings(frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
                              family = wiener(contaminant = TRUE,
                                              max_ndt = mx),
                              data = d))
  ci <- confint(fit)["lambda_(Intercept)", ]
  res <- rbind(res, data.frame(seed = seed, share = mean(hit),
                               est = ci[["est"]], lwr = ci[["lwr"]],
                               upr = ci[["upr"]], code = fit$opt$convergence))
}
res$se <- (res$upr - res$lwr) / 3.92
res$hit <- res$lwr < qlogis(0.05) & qlogis(0.05) < res$upr
cat(sprintf("n %d, contaminants on [%g, %g], max_ndt %s, %d replicates\n",
            n, lo, hi, a[[5]], nrep))
cat(sprintf("logit lambda: mean %.4f (truth %.4f), sd %.4f, mean se %.4f\n",
            mean(res$est), qlogis(0.05), sd(res$est), mean(res$se)))
cat(sprintf("sd / se %.3f; covered %d of %d; code 0 on %d\n",
            sd(res$est) / mean(res$se), sum(res$hit), nrep,
            sum(res$code == 0)))
cat(sprintf("estimate minus realized share, logit: mean %.4f, sd %.4f\n",
            mean(res$est - qlogis(res$share)),
            sd(res$est - qlogis(res$share))))

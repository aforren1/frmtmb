# What the 30 lane-only nonzero convergence codes actually are.
d <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-log/"
a <- read.csv(paste0(d, "falsealarm-base.csv"))
b <- read.csv(paste0(d, "falsealarm-fixed.csv"))
i <- which(a$conv == 0 & b$conv != 0)
cat("lane-only nonzero codes:", length(i), "\n\n")
cat(sprintf("%-9s %5s %12s %12s %10s %12s %12s\n", "arm", "seed",
            "|alpha| base", "|alpha| lane", "d logLik", "ll base",
            "ll lane"))
for (k in i) {
  cat(sprintf("%-9s %5d %12.3g %12.3g %10.4f %12.4f %12.4f\n",
              a$arm[k], a$seed[k], abs(a$alpha[k]), abs(b$alpha[k]),
              b$ll[k] - a$ll[k], a$ll[k], b$ll[k]))
}
cat("\nof the", length(i), "lane-only codes:\n")
cat("  |alpha| above 1e3 on the lane:", sum(abs(b$alpha[i]) > 1e3), "\n")
cat("  lane logLik HIGHER than base:", sum(b$ll[i] > a$ll[i] + 1e-8), "\n")
cat("  lane logLik LOWER than base:", sum(b$ll[i] < a$ll[i] - 1e-8), "\n")
cat("  median gain:", format(stats::median(b$ll[i] - a$ll[i]), digits = 6),
    "\n")
cat("\nall 200 fits, |alpha| > 1e3 against a nonzero code, LANE:\n")
print(table(runaway = abs(b$alpha) > 1e3, conv_nonzero = b$conv != 0))
cat("\nall 200 fits, same table on BASE:\n")
print(table(runaway = abs(a$alpha) > 1e3, conv_nonzero = a$conv != 0))

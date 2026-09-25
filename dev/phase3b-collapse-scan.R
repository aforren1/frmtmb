# Item 3.5: how often does lambda collapse to its edge on CLEAN data,
# one subject of 600 trials, mu 0.8, bs 1.4, ndt 0.3, bias 0.5? Seeds
# 1 to 40. Output: dev/phase3b-log/collapse-scan.txt
.libPaths(c(Sys.getenv("P3B_LIB", "C:/Users/adf44/source/r/phase3b-lib"),
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
rows <- list()
for (seed in 1:40) {
  set.seed(seed)
  d <- ddm_simulate(600, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  f <- bf(rt | dec(upper) ~ 1, bias = 0.5)
  fp <- frm(f, family = wiener(), data = d)
  fc <- suppressWarnings(frm(f, family = wiener(contaminant = TRUE),
                             data = d))
  eta <- fixef_by_dpar(fc)$lambda[[1]]
  dg <- paste(utils::capture.output(diagnose(fc)), collapse = " ")
  rows[[seed]] <- data.frame(
    seed = seed, eta = eta, lambda = plogis(eta),
    dll = as.numeric(logLik(fc)) - as.numeric(logLik(fp)),
    flagged = grepl("end of its link: lambda", dg))
}
x <- do.call(rbind, rows)
sink("dev/phase3b-log/collapse-scan.txt")
print(x, digits = 4)
cat("\ncollapsed (eta < -8):", sum(x$eta < -8), "of", nrow(x), "\n")
cat("flagged by diagnose():", sum(x$flagged), "of", nrow(x), "\n")
cat("flagged among collapsed:", sum(x$flagged & x$eta < -8), "\n")
cat("min logLik gain over plain:", format(min(x$dll)), "\n")
sink()
cat(readLines("dev/phase3b-log/collapse-scan.txt"), sep = "\n")

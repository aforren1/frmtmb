# Summarize dev/moparity-rp-gate.R output: per file, the hyperparameter
# rows with |z| > 2.5 and the ESS distribution per variable.
#   Rscript dev/moparity-rp-gate-sum.R <glob>
args <- commandArgs(TRUE)
fs <- Sys.glob(args[1])
z <- do.call(rbind, lapply(fs, function(f) {
  d <- readRDS(f)
  d$file <- basename(f)
  d
}))
z$hyper <- !grepl("^r_", z$var)
cat("files:", length(fs), " rows:", nrow(z), "\n")
h <- z[z$hyper, ]
agg <- aggregate(cbind(ess_c, ess_n, az_mean = abs(z_mean),
                       az_sd = abs(z_sd)) ~ var + iter, data = h,
                 FUN = function(v) round(stats::median(v), 2))
cat("median per hyperparameter:\n")
print(agg)
mx <- aggregate(cbind(az_mean = abs(z_mean), az_sd = abs(z_sd)) ~ var + iter,
                data = h, FUN = function(v) round(max(v), 2))
cat("max per hyperparameter:\n")
print(mx)
cat("rows with |z| > 3:\n")
big <- z[abs(z$z_mean) > 3 | abs(z$z_sd) > 3, ]
print(big[, c("file", "seed", "iter", "var", "z_mean", "z_sd", "ess_c",
              "ess_n")], digits = 3)
per <- aggregate(cbind(mz = pmax(abs(z_mean), abs(z_sd))) ~ file + seed + iter,
                 data = z, FUN = max)
cat("max |z| over all stats per run, quantiles by iter:\n")
for (it in unique(per$iter)) {
  cat(it, ":", round(stats::quantile(per$mz[per$iter == it],
                                     c(0, 0.5, 0.9, 1)), 2), "\n")
}
per$arm <- ifelse(grepl("-base-", per$file), "base", "lane")
cat("per arm: runs, max |z|, runs with max |z| >= 5\n")
for (a in unique(per$arm)) {
  cat(a, sum(per$arm == a), round(max(per$mz[per$arm == a]), 2),
      sum(per$mz[per$arm == a] >= 5), "\n")
}
cat("min ESS (centered, non-centered):", round(min(z$ess_c), 1),
    round(min(z$ess_n), 1), "\n")
u <- 1 / sqrt(1 / z$ess_c + 1 / z$ess_n)
cat("z units per posterior sd of a mean shift, hyperparameters: min",
    round(min(u[z$hyper]), 1), "median", round(stats::median(u[z$hyper]), 1),
    "\n")

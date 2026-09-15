# Lane `learnhier`: how much does the fit KNOW about each variance
# component? Arithmetic on the records, no fits.
#
# The non-Gaussianity explanation for the `sd(ndt)` bias was measured
# and refuted (dev/learnhier-gauss.R). The next candidate is that the
# data carry much less information about the `ndt` deviations than about
# the others, which is what attenuates a maximum likelihood variance
# component. `confint_varcorr()` reports a Wald interval on the LOG
# scale for a standard deviation, so the interval's own width recovers
# that standard error without refitting anything:
#
#   se_log = (log(upr) - log(lwr)) / (2 * qnorm(0.975))
#
# If `sd(ndt)`'s relative standard error is in line with the others,
# limited information is not the explanation either and the cause is
# still open.
source("dev/learnhier-env.R")
fs <- list.files("dev/learnhier-rec", pattern = "^rlddm-[0-9]+[.]rds$",
                 full.names = TRUE)
R <- lapply(fs, readRDS)
z <- stats::qnorm(0.975)
nm <- c("sd_alpha", "sd_drift", "sd_bs", "sd_ndt")
out <- do.call(rbind, lapply(nm, function(k) {
  e <- vapply(R, function(r) r$vc[[paste0(k, ".est")]], numeric(1))
  lo <- vapply(R, function(r) r$vc[[paste0(k, ".lwr")]], numeric(1))
  hi <- vapply(R, function(r) r$vc[[paste0(k, ".upr")]], numeric(1))
  tr <- vapply(R, function(r) r$truth_fitted[[k]], numeric(1))
  se_log <- (log(hi) - log(lo)) / (2 * z)
  data.frame(component = k, n = length(e),
             mean_est = mean(e), mean_realized = mean(tr),
             rel_bias_pct = 100 * mean(e - tr) / mean(tr),
             mean_se_log = mean(se_log),
             stringsAsFactors = FALSE)
}))
out[-(1:2)] <- lapply(out[-(1:2)], function(x) round(x, 4))
print(out, row.names = FALSE)

cat("\ndo the two rlddm failures move together across replicates?\n")
e_sd <- vapply(R, function(r) r$vc[["sd_ndt.est"]], numeric(1))
t_sd <- vapply(R, function(r) r$truth_fitted[["sd_ndt"]], numeric(1))
e_cr <- vapply(R, function(r) r$vc[["cor_bs~ndt.est"]], numeric(1))
t_cr <- vapply(R, function(r) r$truth_fitted[["cor_bs~ndt"]], numeric(1))
cat("  correlation of the two errors over ", length(R),
    " replicates: ", round(stats::cor(e_sd - t_sd, e_cr - t_cr), 3),
    "\n", sep = "")
cat("  sd(ndt) error: mean ", round(mean(e_sd - t_sd), 4), ", never",
    " positive on ", sum(e_sd - t_sd < 0), " of ", length(R), "\n",
    sep = "")
cat("  cor(bs,ndt) error: mean ", round(mean(e_cr - t_cr), 4),
    ", never positive on ", sum(e_cr - t_cr < 0), " of ", length(R),
    "\n", sep = "")

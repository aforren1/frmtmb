# Checking the INSTRUMENT before believing the measurement.
#
#   Rscript dev/rlddm-scripts/rlddm-popndt.R
#
# The learn-rlddm scale row's population non-decision time, reported as
# `plogis(b0) * mean(floors)`, comes back 9.3 ms above the truth the
# design states, which is 5.3 standard errors. This decomposes that gap
# from the row's own recorded fields (dev/rlddm-scripts/
# rlddm-scale-after.tsv) and asks how much of it is the fit.
#
# No seed: it is arithmetic on numbers the row already produced.

r <- list(frac = 0.828719, pop = 0.25929, se = 0.00175853,
          sub_mean = 0.257234, sub_true = 0.254744,
          sd_nat = 0.0333713, sd_true = 0.0393041,
          rmse_ms = 13.9959, cor = 0.940379)
truth <- 0.25
sd_log <- 0.15

ms <- function(x) sprintf("%+.2f ms", 1000 * x)
cat(sprintf("mean(floors), implied by pop/frac : %.6f\n",
            r$pop / r$frac))
cat("\nthe gap, step by step:\n")
edraw <- truth * exp(sd_log^2 / 2)
cat(sprintf("  stated truth                      %.6f\n", truth))
cat(sprintf("  E[ndt] = 0.25 * exp(sd^2/2)       %.6f  %s\n",
            edraw, ms(edraw - truth)))
cat(sprintf("  realized mean of the 100 draws    %.6f  %s\n",
            r$sub_true, ms(r$sub_true - edraw)))
cat(sprintf("  fitted mean of the 100 learners   %.6f  %s\n",
            r$sub_mean, ms(r$sub_mean - r$sub_true)))
cat(sprintf("  plogis(b0) * mean(floors)         %.6f  %s\n",
            r$pop, ms(r$pop - r$sub_mean)))
cat(sprintf("\n  z of the last against the first   %.2f\n",
            abs(r$pop - truth) / r$se))
cat(sprintf("  the fit's own error, as a fraction of the truth's own\n"))
cat(sprintf("  spread                            %.4f\n",
            abs(r$sub_mean - r$sub_true) / r$sd_true))
cat(sprintf("  per-learner RMSE / sd(truth)      %.4f\n",
            (r$rmse_ms / 1000) / r$sd_true))
cat(sprintf("  correlation with the truth        %.4f\n", r$cor))

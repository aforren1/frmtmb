# Reviewer 2, item 6: the random-effect draws the worker's harness makes
# (set.seed(seed); u <- rnorm(30, 0, 0.35), dev/phase3b-eam-recovery7.R)
# are the SAME in every arm for a given seed. How dispersed is mean(u)
# over seeds 1..80, and what would an IDEAL estimator's interval,
# 0.4 + mean(u) +- 1.96 * sd(u) / sqrt(30), cover there?
ideal <- function(seeds) {
  m <- s <- numeric(length(seeds))
  for (i in seq_along(seeds)) {
    set.seed(seeds[i]); u <- rnorm(30, 0, 0.35); m[i] <- mean(u); s[i] <- sd(u)
  }
  se_ml <- s * sqrt(29 / 30) / sqrt(30)
  c(n = length(seeds), sd_mean_u = sd(m), ratio_to_theory = sd(m) / (0.35 / sqrt(30)),
    cover_known_sd = mean(abs(m) < 1.96 * 0.35 / sqrt(30)),
    cover_ml_sd = mean(abs(m) < 1.96 * se_ml))
}
options(width = 150, digits = 4)
print(rbind(`1..24` = ideal(1:24), `1..40` = ideal(1:40), `1..60` = ideal(1:60),
            `1..80` = ideal(1:80), `1..100` = ideal(1:100), `1001..2000` = ideal(1001:2000),
            `1..10000` = ideal(1:10000)))
# the p-value of the 1..60 dispersion under the model: var ratio ~ chi2(59)/59
m <- vapply(1:60, function(s) { set.seed(s); mean(rnorm(30, 0, 0.35)) }, 0)
cat(sprintf("seeds 1..60: sum(m^2)/(0.35^2/30) = %.1f on 60 df, upper p = %.4f\n",
            sum(m^2) / (0.35^2 / 30), pchisq(sum(m^2) / (0.35^2 / 30), 60, lower.tail = FALSE)))
m <- vapply(1:80, function(s) { set.seed(s); mean(rnorm(30, 0, 0.35)) }, 0)
cat(sprintf("seeds 1..80: sum(m^2)/(0.35^2/30) = %.1f on 80 df, upper p = %.4f\n",
            sum(m^2) / (0.35^2 / 30), pchisq(sum(m^2) / (0.35^2 / 30), 80, lower.tail = FALSE)))

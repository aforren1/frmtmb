# Reviewer 2, item 6: on seeds 1..26, which seeds miss? Base plain fit
# (this review), the oracle, and the worker's lane arms on the SAME seeds.
dirs <- c("dev/phase3b-log/recov5", "dev/phase3b-log/recov4", "dev/phase3b-log/recov3")
arm_ci <- function(arm, s) {
  for (dd in dirs) {
    f <- file.path(dd, sprintf("%s-%d.rds", arm, s)); if (!file.exists(f)) next
    r <- readRDS(f); if (!is.null(r$error)) next
    if (!is.matrix(r$confint)) return(NA)
    ci <- r$confint["(Intercept)", ]; return(ci[["lwr"]] < 0.4 && 0.4 < ci[["upr"]])
  }
  NA
}
fs <- list.files("dev/phase3b-review2/cov-base", "^seed-[0-9]+[.]rds$", full.names = TRUE)
seeds <- sort(as.integer(gsub("[^0-9]", "", basename(fs))))
seeds <- seeds[seeds < 1000]
tab <- t(sapply(seeds, function(s) {
  r <- readRDS(sprintf("dev/phase3b-review2/cov-base/seed-%d.rds", s))
  ci <- r$ci["(Intercept)", ]
  set.seed(s); u <- rnorm(30, 0, 0.35)
  c(seed = s, base = ci[["lwr"]] < 0.4 && 0.4 < ci[["upr"]],
    oracle = abs(mean(u)) < 1.959964 * sd(u) * sqrt(29 / 30) / sqrt(30),
    cleft = arm_ci("cleft", s), cens = arm_ci("cens", s), contfix = arm_ci("contfix", s),
    contdl = arm_ci("contdl", s))
}))
print(tab[tab[, "base"] == 0 | tab[, "oracle"] == 0, ])
cat("base covers", sum(tab[, "base"]), "of", nrow(tab), "; oracle", sum(tab[, "oracle"]), "\n")

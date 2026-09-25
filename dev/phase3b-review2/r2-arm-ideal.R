# Reviewer 2, item 6: for each of the worker's arms, the drift
# intercept's Wald coverage on the fits the worker counted
# (dev/phase3b-summarise2.R's record selection), beside what an ORACLE
# interval covers on the SAME seeds: 0.4 + mean(u) +- 1.96 * SE_ml,
# SE_ml = sd(u) sqrt(29/30) / sqrt(30), from the u the harness drew
# (set.seed(seed); u <- rnorm(30, 0, 0.35)). The oracle knows every
# subject's drift exactly, so it has no estimation error at all.
dirs <- c("dev/phase3b-log/recov5", "dev/phase3b-log/recov4", "dev/phase3b-log/recov3")
oracle <- function(seed) {
  set.seed(seed); u <- rnorm(30, 0, 0.35)
  c(m = mean(u), se = sd(u) * sqrt(29 / 30) / sqrt(30))
}
all_k <- all_n <- all_o <- 0
seen <- integer(0)
for (arm in c("cleft", "cens", "contfix", "contdl", "cont", "collapse")) {
  files <- unlist(lapply(dirs, list.files, pattern = sprintf("^%s-[0-9]+[.]rds$", arm),
                         full.names = TRUE))
  seeds <- sort(unique(as.integer(sub(".*-([0-9]+)[.]rds$", "\\1", files))))
  seeds <- seeds[seeds < 9000]
  k <- o <- n <- 0; used <- integer(0)
  for (s in seeds) {
    got <- NULL
    for (dd in dirs) {
      f <- file.path(dd, sprintf("%s-%d.rds", arm, s))
      if (!file.exists(f)) next
      r <- readRDS(f)
      if (!is.null(r$error)) next
      got <- r; break
    }
    if (is.null(got) || !is.matrix(got$confint) ||
        !"(Intercept)" %in% rownames(got$confint)) next
    ci <- got$confint["(Intercept)", ]
    if (!all(is.finite(ci))) next
    n <- n + 1; used <- c(used, s)
    k <- k + (ci[["lwr"]] < 0.4 && 0.4 < ci[["upr"]])
    oc <- oracle(s); o <- o + (abs(oc[["m"]]) < 1.959964 * oc[["se"]])
  }
  all_k <- all_k + k; all_n <- all_n + n; all_o <- all_o + o; seen <- c(seen, used)
  cat(sprintf("%-9s fits %3d  Wald covers %3d (%5.1f%%)  ORACLE on the same seeds %3d (%5.1f%%)  seeds %d..%d\n",
              arm, n, k, 100 * k / n, o, 100 * o / n, min(used), max(used)))
}
cat(sprintf("pooled    fits %3d  Wald covers %3d (%5.1f%%)  ORACLE %3d (%5.1f%%)\n",
            all_n, all_k, 100 * all_k / all_n, all_o, 100 * all_o / all_n))
cat(sprintf("distinct seeds behind the %d pooled fits: %d\n", all_n, length(unique(seen))))

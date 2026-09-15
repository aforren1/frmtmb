res <- readRDS("dev/sgrev-out/loadcost.rds")
cat("per-arm MINIMA over", nrow(res), "rounds (robust to machine load)\n")
m <- apply(res, 2, min)
for (a in names(m)) cat(sprintf("  %-24s %8.4f\n", a, m[a]))
cat("\ndifferences on the minima:\n")
p <- list(c("sample_FIX","sample_BASE"), c("core_FIX","core_BASE"),
          c("POSCTRL_BASE_plus_10ms","sample_BASE"))
for (q in p) cat(sprintf("  %-38s %+8.4f s\n",
  paste(q[1],"-",q[2]), m[q[1]] - m[q[2]]))
cat("\npermutation test on the paired minima, B = 20000, seed 20260915\n")
set.seed(20260915)
for (q in p) {
  d <- res[, q[1]] - res[, q[2]]
  obs <- min(res[, q[1]]) - min(res[, q[2]])
  s <- replicate(20000, {
    f <- sample(c(TRUE, FALSE), nrow(res), TRUE)
    a <- ifelse(f, res[, q[1]], res[, q[2]])
    b <- ifelse(f, res[, q[2]], res[, q[1]])
    min(a) - min(b)
  })
  cat(sprintf("  %-38s obs %+8.4f  p = %.4f\n", paste(q[1],"-",q[2]),
              obs, mean(abs(s) >= abs(obs))))
}

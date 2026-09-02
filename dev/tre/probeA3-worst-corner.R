# Probe A3. Probe B turned up a Laplace log-likelihood error ten times
# the size of anything in A1, at n = 2 and a latent scale well below the
# residual SD. A1 varied nu and the group size; the missing axis is the
# RATIO s/sigma, which sets how much of the per-group conditional's
# shape the latent density itself contributes.
#
# The error is measured two ways, because only one of them matters:
#   * at a fixed parameter vector (the size of the approximation), and
#   * as the displacement of the ARGMAX (what an ML user actually pays).
#
# Run: Rscript dev/tre/probeA3-worst-corner.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeA3.txt")
t0 <- Sys.time()

cat("Probe A3: the Laplace error over the s/sigma ratio\n")
cat("==================================================\n\n")

cat("A. logLik error at the true parameters, 40 groups, sigma = 1.\n")
cat("   Per-group figures, so the columns are comparable.\n\n")
cat(sprintf("%-6s %-4s %s\n", "nu", "n",
            paste(sprintf("%12s", paste0("s=", c(0.1, 0.25, 0.5, 1, 2, 4))),
                  collapse = "")))
grid <- NULL
for (nu in c(2.5, 3, 5, 10)) {
  for (n in c(2, 5, 15)) {
    row <- numeric(0)
    for (s in c(0.1, 0.25, 0.5, 1, 2, 4)) {
      d <- sim_tre(G = 40, n = n, s = s, nu = nu, seed = 21)
      st <- suff_stats(d, c(1, 0.5))
      ex <- aghq_loglik(d, c(1, 0.5), 1, s, nu, K = 101, st = st)
      la <- laplace_loglik(d, c(1, 0.5), 1, s, nu, st = st)
      row <- c(row, (la - ex) / 40)
      grid <- rbind(grid, data.frame(nu = nu, n = n, s = s,
                                     err = la - ex))
    }
    cat(sprintf("%-6s %-4d %s\n", nu, n,
                paste(sprintf("%12.5f", row), collapse = "")))
  }
}
saveRDS(grid, "dev/tre/probeA3-grid.rds")

cat("\n\nB. What that costs the ESTIMATES. Worst corner only:\n")
cat("   n = 2, s = 0.25 and s = 0.1, 100 replicates, exact ML by\n")
cat("   AGHQ against Laplace ML on identical data.\n\n")
cat(sprintf("%-6s %-6s %16s %16s %16s %10s\n", "nu", "s",
            "d sigma", "d scale s", "d beta1", "d logLik"))
worst <- NULL
for (nu in c(2.5, 3, 5, 10)) {
  for (s in c(0.1, 0.25)) {
    acc <- NULL
    for (r in 1:100) {
      d <- sim_tre(G = 40, n = 2, s = s, nu = nu, seed = 7000 + r)
      fl <- fit_laplace(d, nu = nu)
      fe <- fit_exact(d, nu = nu, K = 51,
                      start = c(fl$par[1], fl$par[2],
                                log(fl$par[3]), log(fl$par[4])))
      acc <- rbind(acc, data.frame(
        d_sigma = fl$par[[3]] - fe$par[[3]],
        d_s = fl$par[[4]] - fe$par[[4]],
        d_b1 = fl$par[[2]] - fe$par[[2]],
        rel_s = (fl$par[[4]] - fe$par[[4]]) / fe$par[[4]],
        ex_s = fe$par[[4]], lap_s = fl$par[[4]],
        d_ll = fl$loglik - fe$loglik))
    }
    worst <- rbind(worst, cbind(nu = nu, s = s, acc))
    f <- function(v) sprintf("%8.5f(%6.5f)", mean(v), sd(v))
    cat(sprintf("%-6s %-6.2f %16s %16s %16s %10.4f\n", nu, s,
                f(acc$d_sigma), f(acc$d_s), f(acc$d_b1),
                mean(acc$d_ll)))
  }
}
saveRDS(worst, "dev/tre/probeA3-worst.rds")

cat("\nC. The same as a relative error in the latent scale, and as a\n")
cat("   fraction of that estimate's own sampling sd.\n\n")
cat(sprintf("%-6s %-6s %14s %14s %14s\n", "nu", "s", "mean rel err",
            "sd(exact s)", "bias / sd"))
for (nu in unique(worst$nu)) for (s in unique(worst$s)) {
  z <- worst[worst$nu == nu & worst$s == s, ]
  cat(sprintf("%-6s %-6.2f %13.2f%% %14.5f %13.1f%%\n", nu, s,
              100 * mean(z$rel_s), sd(z$ex_s),
              100 * mean(z$d_s) / sd(z$ex_s)))
}

cat("\nelapsed: ", format(Sys.time() - t0), "\n")
sink()
cat("done\n")

# Probe A2. The number that decides the feature: how far the Laplace ML
# estimates sit from the EXACT (adaptive-quadrature) ML estimates of the
# same model, as a function of nu and of the group size.
#
# A1 showed the log-likelihood error. A likelihood offset that is nearly
# constant in the parameters costs nothing; what matters is whether the
# argmax moves. Both fits maximize the same target here, one
# approximately and one exactly, so the difference isolates the Laplace
# approximation from sampling error, and each replicate is a paired
# comparison on identical data.
#
# Run: Rscript dev/tre/probeA2-parameter-bias.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeA2.txt")
t0 <- Sys.time()

NREP <- 100L
res <- NULL
for (nu in c(2.5, 3, 5, 10, 30)) {
  for (n in c(3, 8, 25)) {
    for (r in seq_len(NREP)) {
      d <- sim_tre(G = 40, n = n, nu = nu, seed = 1000 * r + n)
      fl <- fit_laplace(d, nu = nu)
      fe <- fit_exact(d, nu = nu, K = 51,
                      start = c(fl$par[1], fl$par[2],
                                log(fl$par[3]), log(fl$par[4])))
      res <- rbind(res, data.frame(
        nu = nu, n = n, rep = r,
        lap_b0 = fl$par[[1]], lap_b1 = fl$par[[2]],
        lap_sigma = fl$par[[3]], lap_s = fl$par[[4]],
        ex_b0 = fe$par[[1]], ex_b1 = fe$par[[2]],
        ex_sigma = fe$par[[3]], ex_s = fe$par[[4]],
        lap_ll = fl$loglik, ex_ll = fe$loglik,
        conv = fl$conv + fe$conv))
    }
  }
}
saveRDS(res, "dev/tre/probeA2.rds")

cat("Probe A2: Laplace ML vs exact (AGHQ) ML, ", NREP,
    " replicates, G = 40 groups\n", sep = "")
cat("Truth: beta0 = 1, beta1 = 0.5, sigma = 1, latent scale s = 1\n")
cat("Non-converged fits: ", sum(res$conv != 0), "\n\n", sep = "")

cat("A. Paired difference Laplace - exact, mean (sd) over replicates\n\n")
cat(sprintf("%-6s %-4s %18s %18s %18s %18s\n",
            "nu", "n", "d beta0", "d beta1", "d sigma", "d scale s"))
for (nu in unique(res$nu)) for (n in unique(res$n)) {
  z <- res[res$nu == nu & res$n == n, ]
  f <- function(a, b) sprintf("%9.5f(%7.5f)", mean(z[[a]] - z[[b]]),
                              sd(z[[a]] - z[[b]]))
  cat(sprintf("%-6s %-4d %18s %18s %18s %18s\n", nu, n,
              f("lap_b0", "ex_b0"), f("lap_b1", "ex_b1"),
              f("lap_sigma", "ex_sigma"), f("lap_s", "ex_s")))
}

cat("\nB. The same difference as a percentage of the exact estimate's\n")
cat("   own Monte-Carlo standard error (how much of a standard error\n")
cat("   the approximation costs)\n\n")
cat(sprintf("%-6s %-4s %10s %10s %10s %10s\n",
            "nu", "n", "beta0", "beta1", "sigma", "scale s"))
for (nu in unique(res$nu)) for (n in unique(res$n)) {
  z <- res[res$nu == nu & res$n == n, ]
  f <- function(a, b) sprintf("%9.2f%%",
                              100 * mean(z[[a]] - z[[b]]) / sd(z[[b]]))
  cat(sprintf("%-6s %-4d %10s %10s %10s %10s\n", nu, n,
              f("lap_b0", "ex_b0"), f("lap_b1", "ex_b1"),
              f("lap_sigma", "ex_sigma"), f("lap_s", "ex_s")))
}

cat("\nC. Bias against the truth, both estimators\n\n")
cat(sprintf("%-6s %-4s %22s %22s\n", "nu", "n",
            "sigma: lap / exact", "scale s: lap / exact"))
for (nu in unique(res$nu)) for (n in unique(res$n)) {
  z <- res[res$nu == nu & res$n == n, ]
  cat(sprintf("%-6s %-4d %10.5f %10.5f %10.5f %10.5f\n", nu, n,
              mean(z$lap_sigma) - 1, mean(z$ex_sigma) - 1,
              mean(z$lap_s) - 1, mean(z$ex_s) - 1))
}

cat("\nD. Log-likelihood gap at each fit's own optimum\n\n")
cat(sprintf("%-6s %-4s %14s %14s\n", "nu", "n", "mean gap", "max gap"))
for (nu in unique(res$nu)) for (n in unique(res$n)) {
  z <- res[res$nu == nu & res$n == n, ]
  cat(sprintf("%-6s %-4d %14.5f %14.5f\n", nu, n,
              mean(z$lap_ll - z$ex_ll), max(abs(z$lap_ll - z$ex_ll))))
}

cat("\nelapsed: ", format(Sys.time() - t0), "\n")
sink()
cat("done\n")

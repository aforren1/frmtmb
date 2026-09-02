# Probe D2. Choosing the fixed nu default on evidence.
#
# Probe C showed nu is not reliably estimable by ML, so the package has
# to ship a default. The two things it trades off:
#   * robustness, which is how little the variance component moves
#     when one group is contaminated (smaller nu is better), and
#   * the cost when nothing is wrong, which is how far a t(nu) fit's
#     variance component sits from a gaussian truth (larger nu is
#     better).
#
# Both are measured on the same contaminated design as probe D, over a
# sweep of nu. The reported latent SD is s * sqrt(nu/(nu-2)).
#
# Run: Rscript dev/tre/probeD2-nu-default.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeD2.txt")
t0 <- Sys.time()

cat("Probe D2: which fixed nu to default to\n")
cat("======================================\n\n")

sim_contaminated <- function(G = 40, n = 8, delta = 5, seed = 1) {
  set.seed(seed)
  b <- rnorm(G)
  b[G] <- b[G] + delta
  sim_tre(G = G, n = n, beta = c(1, 0.5), sigma = 1, s = 1, nu = Inf,
          seed = seed, b_override = b)
}
t_sd <- function(s, nu) s * sqrt(nu / (nu - 2))

NREP <- 80L
nus <- c(2.5, 3, 4, 5, 7, 10, 20)
acc <- NULL
for (delta in c(0, 3, 6, 12)) {
  for (r in seq_len(NREP)) {
    d <- sim_contaminated(delta = delta, seed = 4000 + r)
    fg <- fit_gaussian(d)
    row <- data.frame(delta = delta, rep = r, nu = Inf,
                      sd = fg$par[[4]], b0 = fg$par[[1]],
                      b1 = fg$par[[2]], sigma = fg$par[[3]])
    for (nu in nus) {
      f <- fit_laplace(d, nu = nu)
      row <- rbind(row, data.frame(delta = delta, rep = r, nu = nu,
                                   sd = t_sd(f$par[[4]], nu),
                                   b0 = f$par[[1]], b1 = f$par[[2]],
                                   sigma = f$par[[3]]))
    }
    acc <- rbind(acc, row)
  }
}
saveRDS(acc, "dev/tre/probeD2.rds")

rmse <- function(x, tv) sqrt(mean((x - tv)^2))
cat("A. Latent SD (truth 1). Mean over ", NREP, " replicates.\n\n",
    sep = "")
cat(sprintf("%-8s", "delta"))
for (nu in c(nus, Inf)) cat(sprintf("%10s", paste0("nu=", nu)))
cat("\n")
for (delta in unique(acc$delta)) {
  cat(sprintf("%-8.0f", delta))
  for (nu in c(nus, Inf)) {
    cat(sprintf("%10.4f", mean(acc$sd[acc$delta == delta &
                                        acc$nu == nu])))
  }
  cat("\n")
}

cat("\nB. RMSE of the latent SD against the truth 1.\n\n")
cat(sprintf("%-8s", "delta"))
for (nu in c(nus, Inf)) cat(sprintf("%10s", paste0("nu=", nu)))
cat("\n")
for (delta in unique(acc$delta)) {
  cat(sprintf("%-8.0f", delta))
  for (nu in c(nus, Inf)) {
    cat(sprintf("%10.4f", rmse(acc$sd[acc$delta == delta &
                                        acc$nu == nu], 1)))
  }
  cat("\n")
}

cat("\nC. RMSE of beta0 against the truth 1.\n\n")
cat(sprintf("%-8s", "delta"))
for (nu in c(nus, Inf)) cat(sprintf("%10s", paste0("nu=", nu)))
cat("\n")
for (delta in unique(acc$delta)) {
  cat(sprintf("%-8.0f", delta))
  for (nu in c(nus, Inf)) {
    cat(sprintf("%10.4f", rmse(acc$b0[acc$delta == delta &
                                        acc$nu == nu], 1)))
  }
  cat("\n")
}

cat("\nD. Worst-case latent-SD RMSE across the four contamination\n")
cat("   levels - the minimax choice of nu.\n\n")
cat(sprintf("%-10s %12s %12s\n", "nu", "worst RMSE", "RMSE at delta=0"))
for (nu in c(nus, Inf)) {
  w <- max(sapply(unique(acc$delta), function(dl)
    rmse(acc$sd[acc$delta == dl & acc$nu == nu], 1)))
  cat(sprintf("%-10s %12.4f %12.4f\n", nu, w,
              rmse(acc$sd[acc$delta == 0 & acc$nu == nu], 1)))
}

cat("\nE. Same, for beta0.\n\n")
cat(sprintf("%-10s %12s %12s\n", "nu", "worst RMSE", "RMSE at delta=0"))
for (nu in c(nus, Inf)) {
  w <- max(sapply(unique(acc$delta), function(dl)
    rmse(acc$b0[acc$delta == dl & acc$nu == nu], 1)))
  cat(sprintf("%-10s %12.4f %12.4f\n", nu, w,
              rmse(acc$b0[acc$delta == 0 & acc$nu == nu], 1)))
}

cat("\nelapsed: ", format(Sys.time() - t0), "\n")
sink()
cat("done\n")

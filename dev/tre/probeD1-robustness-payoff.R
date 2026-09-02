# Probe D. What the feature is FOR: does a t latent actually protect the
# variance component and the fixed effects from an outlying group?
#
# The design is the robustlmm motivation stated as a measurement. `G-1`
# groups are drawn from a clean gaussian latent with SD 1; one further
# group is displaced by `delta`. The gaussian-latent fit has to inflate
# the variance component to cover that group, and the fixed effects
# follow it. The t-latent fit should be able to leave it in the tail.
#
# Run: Rscript dev/tre/probeD1-robustness-payoff.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeD1.txt")
t0 <- Sys.time()

cat("Probe D: the robustness payoff, measured\n")
cat("========================================\n\n")
cat("Truth: 39 clean groups with b ~ N(0, 1), one contaminated group\n")
cat("shifted by `delta`; n = 8 per group, sigma = 1, beta = (1, 0.5).\n")
cat("The t fit reports the latent SCALE s; the SD it implies is\n")
cat("s * sqrt(nu / (nu - 2)), which is what compares with the\n")
cat("gaussian fit's SD.\n\n")

sim_contaminated <- function(G = 40, n = 8, delta = 5, seed = 1) {
  set.seed(seed)
  b <- rnorm(G)
  b[G] <- b[G] + delta
  sim_tre(G = G, n = n, beta = c(1, 0.5), sigma = 1, s = 1, nu = Inf,
          seed = seed, b_override = b)
}

t_sd <- function(s, nu) s * sqrt(nu / (nu - 2))

NREP <- 60L
res <- NULL
for (delta in c(0, 2, 5, 10)) {
  for (r in seq_len(NREP)) {
    d <- sim_contaminated(delta = delta, seed = 4000 + r)
    fg <- fit_gaussian(d)
    ft3 <- fit_laplace(d, nu = 3)
    ft5 <- fit_laplace(d, nu = 5)
    ftf <- tryCatch(fit_laplace(d, nu = 5, estimate_nu = TRUE),
                    error = function(e) NULL)
    res <- rbind(res, data.frame(
      delta = delta, rep = r,
      g_b0 = fg$par[[1]], g_b1 = fg$par[[2]],
      g_sigma = fg$par[[3]], g_sd = fg$par[[4]],
      t3_b0 = ft3$par[[1]], t3_b1 = ft3$par[[2]],
      t3_sigma = ft3$par[[3]], t3_s = ft3$par[[4]],
      t3_sd = t_sd(ft3$par[[4]], 3),
      t5_b0 = ft5$par[[1]], t5_b1 = ft5$par[[2]],
      t5_sigma = ft5$par[[3]], t5_s = ft5$par[[4]],
      t5_sd = t_sd(ft5$par[[4]], 5),
      tf_nu = if (is.null(ftf)) NA else ftf$par[[5]],
      tf_sd = if (is.null(ftf)) NA else
        t_sd(ftf$par[[4]], max(ftf$par[[5]], 2.001))))
  }
}
saveRDS(res, "dev/tre/probeD1.rds")

cat("A. Latent SD of the CLEAN process (truth 1) as the contaminating\n")
cat("   group moves away. Mean over ", NREP, " replicates.\n\n", sep = "")
cat(sprintf("%-8s %10s %10s %10s %14s %10s\n", "delta", "gaussian",
            "t(nu=3)", "t(nu=5)", "t(nu free)", "nu_hat"))
for (delta in unique(res$delta)) {
  z <- res[res$delta == delta, ]
  cat(sprintf("%-8.1f %10.4f %10.4f %10.4f %14.4f %10.2f\n", delta,
              mean(z$g_sd), mean(z$t3_sd), mean(z$t5_sd),
              mean(z$tf_sd, na.rm = TRUE), median(z$tf_nu, na.rm = TRUE)))
}

cat("\nB. Latent SCALE the t fits report (what VarCorr would print if\n")
cat("   the scale were reported unconverted).\n\n")
cat(sprintf("%-8s %10s %10s\n", "delta", "t(3) s", "t(5) s"))
for (delta in unique(res$delta)) {
  z <- res[res$delta == delta, ]
  cat(sprintf("%-8.1f %10.4f %10.4f\n", delta, mean(z$t3_s),
              mean(z$t5_s)))
}

cat("\nC. Residual SD (truth 1).\n\n")
cat(sprintf("%-8s %10s %10s %10s\n", "delta", "gaussian", "t(3)", "t(5)"))
for (delta in unique(res$delta)) {
  z <- res[res$delta == delta, ]
  cat(sprintf("%-8.1f %10.4f %10.4f %10.4f\n", delta, mean(z$g_sigma),
              mean(z$t3_sigma), mean(z$t5_sigma)))
}

cat("\nD. Fixed effects: RMSE against the truth (beta0 = 1,",
    "beta1 = 0.5).\n\n")
cat(sprintf("%-8s %-8s %10s %10s %10s\n", "delta", "coef", "gaussian",
            "t(3)", "t(5)"))
rmse <- function(x, tv) sqrt(mean((x - tv)^2))
for (delta in unique(res$delta)) {
  z <- res[res$delta == delta, ]
  cat(sprintf("%-8.1f %-8s %10.4f %10.4f %10.4f\n", delta, "beta0",
              rmse(z$g_b0, 1), rmse(z$t3_b0, 1), rmse(z$t5_b0, 1)))
  cat(sprintf("%-8.1f %-8s %10.4f %10.4f %10.4f\n", delta, "beta1",
              rmse(z$g_b1, 0.5), rmse(z$t3_b1, 0.5), rmse(z$t5_b1, 0.5)))
}

cat("\nE. The cost when there is NO contamination (delta = 0): the t\n")
cat("   fit's efficiency loss relative to the correctly specified\n")
cat("   gaussian fit, as a ratio of RMSEs (>1 means the t costs).\n\n")
z <- res[res$delta == 0, ]
cat(sprintf("beta0   t(3)/gaussian %.4f   t(5)/gaussian %.4f\n",
            rmse(z$t3_b0, 1) / rmse(z$g_b0, 1),
            rmse(z$t5_b0, 1) / rmse(z$g_b0, 1)))
cat(sprintf("beta1   t(3)/gaussian %.4f   t(5)/gaussian %.4f\n",
            rmse(z$t3_b1, 0.5) / rmse(z$g_b1, 0.5),
            rmse(z$t5_b1, 0.5) / rmse(z$g_b1, 0.5)))
cat(sprintf("sd      t(3)/gaussian %.4f   t(5)/gaussian %.4f\n",
            rmse(z$t3_sd, 1) / rmse(z$g_sd, 1),
            rmse(z$t5_sd, 1) / rmse(z$g_sd, 1)))

cat("\nF. One dataset in detail (delta = 10, seed 4001): where the\n")
cat("   contaminated group's own predicted latent lands.\n\n")
d <- sim_contaminated(delta = 10, seed = 4001)
fg <- fit_gaussian(d)
ft <- fit_laplace(d, nu = 3)
cat(sprintf("true b[40]            %8.4f\n", d$b[40]))
cat(sprintf("gaussian b_hat[40]    %8.4f   latent SD %.4f\n",
            fg$b[40], fg$par[[4]]))
cat(sprintf("t(3)     b_hat[40]    %8.4f   latent SD %.4f (scale %.4f)\n",
            ft$b[40], t_sd(ft$par[[4]], 3), ft$par[[4]]))
cat(sprintf("gaussian b_hat[1:5]   %s\n",
            paste(sprintf("%7.3f", fg$b[1:5]), collapse = "")))
cat(sprintf("t(3)     b_hat[1:5]   %s\n",
            paste(sprintf("%7.3f", ft$b[1:5]), collapse = "")))
cat(sprintf("true     b[1:5]       %s\n",
            paste(sprintf("%7.3f", d$b[1:5]), collapse = "")))
cat(sprintf("\nshrinkage of the clean groups, RMSE(b_hat - b_true):\n"))
cat(sprintf("  gaussian %.4f   t(3) %.4f\n",
            sqrt(mean((fg$b[1:39] - d$b[1:39])^2)),
            sqrt(mean((ft$b[1:39] - d$b[1:39])^2))))

cat("\nelapsed: ", format(Sys.time() - t0), "\n")
sink()
cat("done\n")

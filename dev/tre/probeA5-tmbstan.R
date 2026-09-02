# Probe A5. A reference that does not go through my own quadrature.
#
# `tmbstan` runs NUTS over the same RTMB object two ways: `laplace =
# FALSE` samples the FULL joint (fixed parameters and latents together,
# no approximation anywhere), and `laplace = TRUE` samples the
# Laplace-marginalized objective. With flat priors on the transformed
# parameters both target the same posterior, so any systematic
# difference between the two is the Laplace error expressed on the
# posterior scale, measured with no quadrature and no code of mine in
# the path.
#
# n = 3 per group is the worst case A2 found at a realistic scale;
# n = 8 is the ordinary one.
#
# Run: Rscript dev/tre/probeA5-tmbstan.R
source("dev/tre/tre-common.R")
suppressMessages(library(tmbstan))

sink("dev/tre/probeA5.txt")
options(mc.cores = 1)

cat("Probe A5: full-joint NUTS against Laplace NUTS on the same tape\n")
cat("===============================================================\n")
cat("tmbstan ", as.character(packageVersion("tmbstan")),
    ", rstan ", as.character(packageVersion("rstan")), "\n\n", sep = "")

pnames <- c("beta[1]", "beta[2]", "log_sigma", "log_s")
for (nu in c(3, 5)) {
  for (n in c(3, 8)) {
    d <- sim_tre(G = 40, n = n, nu = nu, seed = 77)
    fl <- fit_laplace(d, nu = nu)
    set.seed(1)
    sj <- suppressWarnings(tmbstan(fl$obj, chains = 4, iter = 4000,
                                   warmup = 1000, laplace = FALSE,
                                   seed = 11, refresh = 0,
                                   control = list(adapt_delta = 0.95)))
    set.seed(1)
    sl <- suppressWarnings(tmbstan(fl$obj, chains = 4, iter = 4000,
                                   warmup = 1000, laplace = TRUE,
                                   seed = 11, refresh = 0,
                                   control = list(adapt_delta = 0.95)))
    mj <- as.matrix(sj)[, pnames, drop = FALSE]
    ml <- as.matrix(sl)[, pnames, drop = FALSE]
    ess <- min(rstan::summary(sj)$summary[pnames, "n_eff"],
               rstan::summary(sl)$summary[pnames, "n_eff"])
    rh <- max(rstan::summary(sj)$summary[pnames, "Rhat"],
              rstan::summary(sl)$summary[pnames, "Rhat"])
    cat(sprintf("\nnu = %s, n = %s per group, 40 groups", nu, n),
        sprintf("   (min n_eff %.0f, max Rhat %.3f)\n", ess, rh))
    cat(sprintf("%-12s %10s %10s %10s %10s %10s %12s\n", "parameter",
                "joint mean", "lap mean", "joint sd", "lap sd",
                "Laplace ML", "d/sd"))
    for (p in pnames) {
      mle <- switch(p, "beta[1]" = fl$par[[1]], "beta[2]" = fl$par[[2]],
                    "log_sigma" = log(fl$par[[3]]),
                    "log_s" = log(fl$par[[4]]))
      cat(sprintf("%-12s %10.5f %10.5f %10.5f %10.5f %10.5f %11.2f%%\n",
                  p, mean(mj[, p]), mean(ml[, p]), sd(mj[, p]),
                  sd(ml[, p]), mle,
                  100 * (mean(ml[, p]) - mean(mj[, p])) / sd(mj[, p])))
    }
    # the marginal posterior of the latent for the most extreme group,
    # where the fat tail is doing the most work
    j <- which.max(abs(d$b))
    bj <- paste0("b[", j, "]")
    if (bj %in% colnames(as.matrix(sj))) {
      bs <- as.matrix(sj)[, bj]
      cat(sprintf("%-12s joint posterior mean %8.4f sd %7.4f; ",
                  bj, mean(bs), sd(bs)))
      cat(sprintf("Laplace mode %8.4f; true %8.4f\n", fl$b[j], d$b[j]))
    }
  }
}
sink()
cat("done\n")

# Probe C. Is nu identifiable by maximum likelihood in a realistic
# design?
#
# brms estimates nu (`df_1`, a real<lower=1>) with a gamma(2, 0.1) prior
# truncated at 1. A prior with mean 20 over a parameter the data barely
# see is doing real work there; frmtmb has no prior to lean on, so the
# question is what the LIKELIHOOD alone says about nu.
#
# C1 profiles the likelihood in nu. C2 estimates nu jointly and counts
# how often it runs to a boundary. C3 checks that RTMB can tape `dt`
# with an AD degrees-of-freedom argument at all. C4 asks the easier
# question that matters more in practice: can the likelihood tell a t
# latent from a gaussian one?
#
# Run: Rscript dev/tre/probeC1-nu-identifiability.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeC1.txt")
t0 <- Sys.time()

cat("Probe C: identifiability of the latent degrees of freedom\n")
cat("=========================================================\n\n")

# --- C3 first: can nu even be a parameter? --------------------------
cat("C3. RTMB tapes dt() with an AD `df`.\n\n")
d <- sim_tre(G = 40, n = 5, nu = 3, seed = 3)
fl <- fit_laplace(d, nu = 3, estimate_nu = TRUE)
g_ad <- fl$obj$gr(fl$opt$par)
eps <- 1e-5
g_fd <- sapply(seq_along(fl$opt$par), function(i) {
  p1 <- p2 <- fl$opt$par
  p1[i] <- p1[i] + eps; p2[i] <- p2[i] - eps
  (fl$obj$fn(p1) - fl$obj$fn(p2)) / (2 * eps)
})
cat("AD gradient  : ", paste(sprintf("%12.3e", g_ad), collapse = ""), "\n")
cat("FD gradient  : ", paste(sprintf("%12.3e", g_fd), collapse = ""), "\n")
cat("max abs diff : ", sprintf("%.3e", max(abs(g_ad - g_fd))), "\n")
cat("nu_hat       : ", sprintf("%.4f", fl$par[["nu"]]), " (true 3)\n\n")

# --- C1: the profile ------------------------------------------------
cat("\nC1. Profile log-likelihood in nu (Laplace), true nu = 3,\n")
cat("    n = 5 observations per group. Reported as the drop from the\n")
cat("    profile maximum; a 1.92 drop is the 95% profile interval.\n\n")
nugrid <- c(2.1, 2.5, 3, 4, 5, 7, 10, 15, 25, 50, 100, 500)
prof <- NULL
for (G in c(20, 40, 100, 400)) {
  d <- sim_tre(G = G, n = 5, nu = 3, seed = 42)
  lls <- sapply(nugrid, function(nu) fit_laplace(d, nu = nu)$loglik)
  best <- which.max(lls)
  ci <- range(nugrid[lls >= max(lls) - 1.92])
  prof <- rbind(prof, data.frame(G = G, nu = nugrid, ll = lls))
  cat(sprintf("G = %-4d peak nu = %-6s  95%% profile interval [%s, %s]\n",
              G, nugrid[best], ci[1], ci[2]))
  cat("  drop: ")
  cat(paste(sprintf("%s:%.2f", nugrid, max(lls) - lls), collapse = "  "))
  cat("\n\n")
}
saveRDS(prof, "dev/tre/probeC1-profile.rds")

# --- C2: joint estimation over replicates ---------------------------
cat("\nC2. Joint ML for nu, 100 replicates. nu = 2 + exp(log_nu), so\n")
cat("    'at ceiling' means the optimizer walked off toward the\n")
cat("    gaussian limit and 'at floor' means it collapsed on nu = 2.\n\n")
NREP <- 100L
c2 <- NULL
for (true_nu in c(3, 5, 10)) {
  for (G in c(20, 40, 100)) {
    nus <- numeric(NREP)
    for (r in seq_len(NREP)) {
      dd <- sim_tre(G = G, n = 5, nu = true_nu, seed = 5000 + r)
      f <- tryCatch(fit_laplace(dd, nu = true_nu, estimate_nu = TRUE),
                    error = function(e) NULL)
      nus[r] <- if (is.null(f)) NA_real_ else f$par[["nu"]]
    }
    c2 <- rbind(c2, data.frame(true_nu = true_nu, G = G,
                               nu = nus, rep = seq_len(NREP)))
    cat(sprintf(paste0("true nu %-4s G %-4d  median %8.2f  ",
                       "IQR [%7.2f, %8.2f]  ceiling(>100) %5.1f%%  ",
                       "floor(<2.1) %5.1f%%  failed %d\n"),
                true_nu, G, median(nus, na.rm = TRUE),
                quantile(nus, 0.25, na.rm = TRUE),
                quantile(nus, 0.75, na.rm = TRUE),
                100 * mean(nus > 100, na.rm = TRUE),
                100 * mean(nus < 2.1, na.rm = TRUE),
                sum(is.na(nus))))
  }
}
saveRDS(c2, "dev/tre/probeC1-joint.rds")

# --- C4: can the likelihood tell t from gaussian? -------------------
cat("\n\nC4. Detecting a t latent at all: 2*(logLik(nu free) -\n")
cat("    logLik(gaussian)). The boundary-free comparison is against\n")
cat("    nu = 500 rather than a separate gaussian fit.\n\n")
cat(sprintf("%-8s %-5s %12s %12s %12s\n", "true nu", "G",
            "median LR", "P(LR>3.84)", "P(LR>0.01)"))
for (true_nu in c(3, 5, 10, Inf)) {
  for (G in c(20, 40, 100)) {
    lr <- numeric(40)
    for (r in seq_len(40)) {
      dd <- sim_tre(G = G, n = 5, nu = true_nu, seed = 9000 + r)
      f1 <- tryCatch(fit_laplace(dd, nu = 3, estimate_nu = TRUE),
                     error = function(e) NULL)
      f0 <- tryCatch(fit_laplace(dd, nu = 500), error = function(e) NULL)
      lr[r] <- if (is.null(f1) || is.null(f0)) NA_real_ else
        2 * (f1$loglik - f0$loglik)
    }
    lr <- pmax(lr, 0)
    cat(sprintf("%-8s %-5d %12.3f %11.1f%% %11.1f%%\n", true_nu, G,
                median(lr, na.rm = TRUE),
                100 * mean(lr > 3.84, na.rm = TRUE),
                100 * mean(lr > 0.01, na.rm = TRUE)))
  }
}

cat("\nelapsed: ", format(Sys.time() - t0), "\n")
sink()
cat("done\n")

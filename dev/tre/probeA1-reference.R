# Probe A1. Certify the exact reference, then measure the Laplace error
# in the log-likelihood at a FIXED parameter vector (the truth), so that
# the number is the approximation error and nothing else.
#
# Run: Rscript dev/tre/probeA1-reference.R
source("dev/tre/tre-common.R")

log_to <- "dev/tre/probeA1.txt"
sink(log_to)

cat("Probe A1: reference certification and Laplace error at the truth\n")
cat("================================================================\n\n")

cat("1. Does AGHQ agree with integrate() over the whole line?\n\n")
cat(sprintf("%-6s %-6s %14s %14s %14s %14s %12s\n",
            "nu", "n", "AGHQ K=21", "AGHQ K=51", "AGHQ K=101",
            "integrate", "|K101-int|"))
cert <- list()
for (nu in c(2.5, 3, 5, 10, 30)) {
  for (n in c(4, 20)) {
    d <- sim_tre(G = 40, n = n, nu = nu, seed = 11)
    st <- suff_stats(d, c(1, 0.5))
    a21 <- aghq_loglik(d, c(1, 0.5), 1, 1, nu, K = 21, st = st)
    a51 <- aghq_loglik(d, c(1, 0.5), 1, 1, nu, K = 51, st = st)
    a101 <- aghq_loglik(d, c(1, 0.5), 1, 1, nu, K = 101, st = st)
    itg <- integrate_loglik(d, c(1, 0.5), 1, 1, nu, st = st)
    cat(sprintf("%-6s %-6d %14.9f %14.9f %14.9f %14.9f %12.2e\n",
                nu, n, a21, a51, a101, itg, abs(a101 - itg)))
    cert[[paste(nu, n)]] <- c(a21 = a21, a51 = a51, a101 = a101, int = itg)
  }
}

cat("\n2. Laplace error in the log-likelihood at the true parameters.\n")
cat("   err = logLik(Laplace) - logLik(exact); 40 groups.\n\n")
cat(sprintf("%-6s %-4s %14s %14s %12s %12s\n",
            "nu", "n", "exact", "Laplace", "abs err", "rel err"))
tabA <- NULL
for (nu in c(2.5, 3, 4, 5, 10, 30, Inf)) {
  for (n in c(3, 5, 10, 25)) {
    d <- sim_tre(G = 40, n = n, nu = nu, seed = 11)
    st <- suff_stats(d, c(1, 0.5))
    ex <- aghq_loglik(d, c(1, 0.5), 1, 1, nu, K = 101, st = st)
    la <- laplace_loglik(d, c(1, 0.5), 1, 1, nu, st = st)
    cat(sprintf("%-6s %-4d %14.6f %14.6f %12.4f %12.2e\n",
                nu, n, ex, la, la - ex, abs((la - ex) / ex)))
    tabA <- rbind(tabA, data.frame(nu = nu, n = n, exact = ex,
                                   laplace = la, err = la - ex))
  }
}

cat("\n3. Does the hand-rolled Laplace match TMB's?  (nu = 3, n = 8)\n\n")
d <- sim_tre(G = 40, n = 8, nu = 3, seed = 11)
fl <- fit_laplace(d, nu = 3)
# TMB's objective at the SAME parameters as the hand-rolled one
p <- c(beta = c(1, 0.5), log_sigma = 0, log_s = 0)
tmb_at_truth <- -fl$obj$fn(p)
hand_at_truth <- laplace_loglik(d, c(1, 0.5), 1, 1, 3)
cat(sprintf("TMB Laplace at truth   %.9f\n", tmb_at_truth))
cat(sprintf("hand Laplace at truth  %.9f\n", hand_at_truth))
cat(sprintf("difference             %.3e\n",
            abs(tmb_at_truth - hand_at_truth)))
cat(sprintf("exact at truth         %.9f\n",
            aghq_loglik(d, c(1, 0.5), 1, 1, 3, K = 101)))

saveRDS(list(cert = cert, tabA = tabA,
             tmb_check = c(tmb = tmb_at_truth, hand = hand_at_truth)),
        "dev/tre/probeA1.rds")
sink()
cat("done\n")

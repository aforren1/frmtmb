# Reviewer, punch round 2: WHY the two x2 slopes under-cover, since the
# answer decides whether the remedy `?lca` names can work.
#
# `confint(method = "profile")` widens an interval when the profile
# log-likelihood is not quadratic in the parameter. If instead the
# reported standard error is simply below the estimator's sampling
# spread, profiling cannot help, because it is not a curvature problem.
# The z-scores say which: a heavy tail with a normal-sized middle points
# at occasional badly identified fits, and a uniformly inflated spread
# points at the standard error itself.
#
#   Rscript dev/rev-latent-ztail.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env3.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
for (f in c("latent-2p4-recheck-oos.tsv", "latent-2p4-recheck.tsv")) {
  r <- utils::read.delim(paste0(D, f))
  cat("\n==== ", f, " ====\n", sep = "")
  cat(sprintf("%-6s %8s %10s %8s %10s %8s\n", "coef", "sd(z)",
              "IQR(z)/1.349", "ratio", "kurtosis", "|z|>1.96"))
  for (nm in c("c2_3", "c3_3", "c4_3", "c3_2")) {
    e <- r[[paste0("err_", nm)]]
    s <- r[[paste0("se_", nm)]]
    z <- e / s
    rob <- stats::IQR(z) / 1.349
    cat(sprintf("%-6s %8.4f %10.4f %8.4f %10.3f %5d/%d\n", nm,
                stats::sd(z), rob, stats::sd(z) / rob,
                mean((z - mean(z))^4) / stats::sd(z)^4,
                sum(abs(z) > 1.96), length(z)))
  }
}
cat("\nsd(z) is 1 when the reported standard error matches the spread.\n")
cat("sd(z) / (IQR(z)/1.349) near 1 means the inflation is in the whole\n")
cat("distribution and not in a tail, so there is no curvature for a\n")
cat("profile interval to find.\n")

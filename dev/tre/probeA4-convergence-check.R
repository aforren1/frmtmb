# Probe A4. Two loose ends from A2 and A3.
#
# (1) A2 reported a non-zero `nlminb` convergence code on most fits.
#     `rel.tol = 1e-12` is far inside the noise floor of an AGHQ
#     objective evaluated by 51-node quadrature, so the code is
#     expected to be a tolerance artefact rather than a failed fit. The
#     check that settles it is the GRADIENT, not the code.
#
# (2) A3 showed the log-likelihood error approaching a per-group
#     CONSTANT as the latent scale shrinks against the residual SD.
#     That constant has a closed form. It is the Laplace error of the
#     scaled t density on its own, which is scale-free:
#
#       c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) + 0.5*log(2/(nu+1))
#
#     A constant offset in the log-likelihood does not move the argmax,
#     which is why the estimates in A2 and A3 survive an error that
#     looks large in the likelihood. It DOES move logLik/AIC/BIC, so it
#     is the number the help page has to name.
#
# Run: Rscript dev/tre/probeA4-convergence-check.R
source("dev/tre/tre-common.R")

sink("dev/tre/probeA4.txt")

cat("Probe A4: convergence codes and the small-scale constant\n")
cat("========================================================\n\n")

cat("1. Gradient norms at the optima A2 flagged as non-converged.\n\n")
cat(sprintf("%-6s %-4s %10s %14s %14s %s\n", "nu", "n", "conv code",
            "max|grad| lap", "max|grad| exact", "nlminb message"))
for (nu in c(2.5, 3, 10)) {
  for (n in c(3, 8)) {
    d <- sim_tre(G = 40, n = n, nu = nu, seed = 1003 + n)
    fl <- fit_laplace(d, nu = nu)
    fe <- fit_exact(d, nu = nu, K = 51,
                    start = c(fl$par[1], fl$par[2], log(fl$par[3]),
                              log(fl$par[4])))
    gl <- max(abs(fl$obj$gr(fl$opt$par)))
    obj_e <- function(p) -aghq_loglik(d, p[1:2], exp(p[3]), exp(p[4]),
                                      nu, K = 51)
    pe <- c(fe$par[[1]], fe$par[[2]], log(fe$par[[3]]),
            log(fe$par[[4]]))
    ge <- max(abs(numDeriv::grad(obj_e, pe)))
    cat(sprintf("%-6s %-4d %10d %14.2e %14.2e %s\n", nu, n,
                fl$conv + fe$conv, gl, ge,
                paste(unique(c(fl$opt$message, fe$opt$message)),
                      collapse = " / ")))
  }
}

cat("\n2. The small-scale limit constant, closed form against the\n")
cat("   measured per-group error at s/sigma = 0.02.\n\n")
cconst <- function(nu) {
  lgamma((nu + 1) / 2) - lgamma(nu / 2) + 0.5 * log(2 / (nu + 1))
}
cat(sprintf("%-8s %16s %16s %12s\n", "nu", "closed form c(nu)",
            "measured / group", "diff"))
for (nu in c(2.5, 3, 4, 5, 7, 10, 20, 30, 100)) {
  s <- 0.02
  d <- sim_tre(G = 40, n = 5, s = s, nu = nu, seed = 31)
  st <- suff_stats(d, c(1, 0.5))
  ex <- aghq_loglik(d, c(1, 0.5), 1, s, nu, K = 101, st = st)
  la <- laplace_loglik(d, c(1, 0.5), 1, s, nu, st = st)
  cat(sprintf("%-8s %16.6f %16.6f %12.2e\n", nu, cconst(nu),
              (la - ex) / 40, abs(cconst(nu) - (la - ex) / 40)))
}

cat("\n3. What that constant does to AIC. Sum over G groups; the two\n")
cat("   models being compared by AIC would each carry their own.\n\n")
cat(sprintf("%-8s %12s %12s %12s %12s\n", "nu", "c(nu)", "G=20",
            "G=100", "G=1000"))
for (nu in c(2.5, 3, 5, 10, 30)) {
  cc <- cconst(nu)
  cat(sprintf("%-8s %12.5f %12.3f %12.3f %12.3f\n", nu, cc, 20 * cc,
              100 * cc, 1000 * cc))
}

cat("\n4. And the gaussian limit: c(nu) -> 0 as nu -> Inf.\n\n")
for (nu in c(1e2, 1e3, 1e4, 1e6)) {
  cat(sprintf("nu = %-10g  c(nu) = %12.3e\n", nu, cconst(nu)))
}

sink()
cat("done\n")

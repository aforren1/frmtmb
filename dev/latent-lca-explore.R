# Lane `latent`, item 2.4: one replicate, to learn the shapes before
# paying for a hundred. Seed 20260909.
#
#   Rscript dev/latent-lca-explore.R

source("dev/latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
lane_env_report()

source("dev/latent-lca-sim.R")

s <- lca_sim(seed = 20260909L)
cat("\nn =", nrow(s$dd), " J =", ncol(s$dd$Y), " K =", nrow(s$base), "\n")
cat("true class shares:", round(prop.table(table(s$cl)), 4), "\n")

t0 <- Sys.time()
fit <- frm(bf(Y ~ x1 + x2), family = lca(K = 4), data = s$dd)
cat("frm seconds:", as.numeric(difftime(Sys.time(), t0, units = "secs")),
    "\n")
cat("logLik:", format(as.numeric(logLik(fit)), digits = 12), "\n")
print(fixef(fit))
V <- vcov(fit)
cat("vcov dim:", dim(V), "\n")
print(rownames(V))
cat("\nextra_names:", fit$frame$extra_names, "\n")
cat("estimates names:", names(fit$estimates), "\n")

pf <- lca_profiles(fit)
cat("\nprofiles length:", length(pf), " dim1:", dim(pf[[1L]]), "\n")
print(attr(pf, "class_sizes"))

suppressMessages(loadNamespace("poLCA"))
dp <- s$dd
dp$Y <- NULL
for (j in seq_len(ncol(s$Y))) dp[[paste0("I", j)]] <- s$Y[, j]
f <- stats::as.formula(paste0("cbind(",
  paste(paste0("I", seq_len(ncol(s$Y))), collapse = ", "), ") ~ x1 + x2"))
set.seed(1L)
t0 <- Sys.time()
pl <- poLCA::poLCA(f, dp, nclass = 4, nrep = 5, verbose = FALSE,
                   maxiter = 20000, tol = 1e-12)
cat("poLCA seconds:", as.numeric(difftime(Sys.time(), t0, units = "secs")),
    "\n")
cat("poLCA llik:", format(pl$llik, digits = 12), "\n")
cat("poLCA coeff dim:", dim(pl$coeff), "\n")
print(pl$coeff)
cat("poLCA probs[[1]]:\n"); print(pl$probs[[1L]])
cat("poLCA P:", pl$P, "\n")
cat("names(pl):", names(pl), "\n")

# Lane seflags, item 1: the spectrum se_tier3() reads on ported fixture
# 1 (data-helpers:7): the flat directions, tau, each parameter's
# projection onto the flat subspace, and the cosine of the five
# predictions' gradients with that subspace in the unit-diagonal
# metric (the one tau bounds). Meaningful where the fixture converges
# (OpenBLAS 0.3.26 with LAPACK).
#   Rscript dev/seflags-dh7tau.R <lib dir>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb)})
cat("lib", find.package("frmtmb"), "| R", R.home(), "\n")
env <- new.env(parent = asNamespace("frmtmb"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/tests/testthat"
for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
ns <- asNamespace("frmtmb")
withr::local_seed(3L)
fit <- suppressWarnings(suppressMessages(env$brms_fixture(1)))
cat("code", fit$opt$convergence, "\n")
sdr <- suppressWarnings(ns$sdr_of(fit))
print(sdr$se_lost)
h <- fit$cache$hessian_fixed
if (is.null(h)) h <- ns$fit_outer_hessian(fit)
H <- h$H; E <- h$E
nm <- ns$outer_par_names(fit)
lostn <- names(sdr$se_lost)
free <- which(!nm %in% lostn[sdr$se_lost != "flat" |
                             nm[match(lostn, nm)] %in% c("theta_3", "theta_4")])
# the free set se_tier3() decomposes: every parameter not removed as a row
# (the two thetas went as rows: edge and empty)
free <- setdiff(seq_along(nm), match(c("theta_3", "theta_4"), nm))
Hf <- H[free, free]; Df <- sqrt(abs(diag(Hf)))
e <- eigen(Hf / outer(Df, Df), symmetric = TRUE)
Es <- E[free, free] / outer(Df, Df)
enorm <- sqrt(sum(Es^2))
nk <- sqrt(colSums((Es %*% e$vectors)^2))
cat("ev:", signif(e$values, 3), "\n")
cat("noise per direction:", signif(nk, 3), "\n")
k <- length(e$values)
gap <- e$values[k - 1] - e$values[k]
tau <- max(sqrt(1e-9), 10 * enorm / gap)
cat("enorm", signif(enorm, 3), "gap", signif(gap, 3), "tau", signif(tau, 3),
    "\n")
v <- e$vectors[, k]
names(v) <- nm[free]
cat("smallest direction (unit diagonal):\n"); print(signif(v, 4))
proj <- abs(v)
cat("0.5 * max proj:", signif(0.5 * max(proj), 3), "\n")
# the predictions' gradients over the outer coefficients: the fixed rows
X <- as.matrix(fit$frame$linpreds[[1]]$X)
bi <- which(names(fit$opt$par) == "beta")
G <- matrix(0, 5, length(nm)); G[, bi] <- X[1:5, ]
Gs <- G[, free] / rep(Df, each = 5)
cosv <- abs(Gs %*% v) / sqrt(rowSums(Gs^2))
cat("prediction cosines with the full direction (unit-diag metric):",
    signif(cosv, 3), "\n")
# the exact alias in the same metric
cn <- colnames(X)
d <- numeric(length(nm))
d[bi[match("(Intercept)", cn)]] <- 0.046394458
d[bi[match("Age", cn)]] <- -0.190041455
d[bi[grep("^s[(]Age", cn)]] <- 1
ds <- d[free] * Df; ds <- ds / sqrt(sum(ds^2))
cat("angle between the measured direction and the exact alias (sin):",
    signif(sqrt(max(0, 1 - sum(ds * v)^2)), 3), "\n")

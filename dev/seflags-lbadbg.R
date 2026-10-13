# Lane seflags: why the LBA fit of test-lba.R:142 lost every SE on the
# first lane build: the spectrum, each small direction's corrected
# curvature and the step test's verdict.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
ns <- asNamespace("frmtmb")
env <- asNamespace("frmtmb.eam")
set.seed(21)
dat <- env$lba_simulate(150, v = c(2.2, 1.4, 0.8), A = 0.5, k = 0.4, ndt = 0.2)
fit <- suppressWarnings(frm(bf(rt | vint(choice) ~ 1), family = lba(3), data = dat))
p <- fit$opt$par
print(p); cat("grad", signif(fit$obj$gr(p), 3), "\n")
H <- fit$obj$he(p); E <- abs(H - t(H)) / 2; H <- (H + t(H)) / 2
D <- sqrt(abs(diag(H))); e <- eigen(H / outer(D, D), symmetric = TRUE)
big <- max(abs(e$values)); cat("ev rel", signif(e$values / big, 3), "\n")
np <- length(p); keep <- rep(TRUE, np)
for (k in which(e$values <= 1e-3 * big)) {
  c0 <- ns$se_ridge_curvature(fit, p, seq_len(np), e, D, k, keep)
  cr <- ns$se_curvature_real(fit, p, seq_len(np), e$vectors[, k] / D, max(c0, 1e-300))
  cat(sprintf("k %d ev %.3g c0 %.3g thr %.3g step-test %s\n", k, e$values[k], c0, 1e-9 * big, cr))
}
k <- np; vk <- e$vectors[, k]
cat("v(k):", signif(vk, 3), "\n")
cur <- function(z, h) { d <- h * vk / D
  sum((as.numeric(fit$obj$gr(z + d)) - as.numeric(fit$obj$gr(z - d))) / D * vk) / (2 * h) }
cat("curv at p:", signif(sapply(c(1e-4, 1e-3, 1e-2, 1e-1, 1), function(h) cur(p, h)), 3), "\n")
f0 <- fit$obj$fn(p)
cat("f along v:", signif(sapply(c(-3, -1, -0.3, 0.3, 1, 3), function(h) fit$obj$fn(p + h * vk / D) - f0), 3), "\n")

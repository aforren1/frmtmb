# Review of lane sparsese: does spHess(random = FALSE) carry the cross
# block H_rf? Compare with a central difference of the joint gradient.
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
ng <- 100L
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
options(frmtmb.joint_root_max = 20)
fit <- suppressMessages(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g),
                            data = d, REML = TRUE))
sdr <- frmtmb:::sdr_of(fit)
obj <- fit$obj
lp <- obj$env$last.par.best
r <- obj$env$random
f <- setdiff(seq_along(lp), r)
H <- as.matrix(Matrix::forceSymmetric(obj$env$spHess(lp, random = FALSE)))
cat("max |H_rf|:", max(abs(H[r, f])), " max |H_ff|:", max(abs(H[f, f])), "\n")
gfull <- function(x) as.numeric(obj$env$f(x, order = 1))
h <- 1e-5
Hfd <- sapply(f, function(k) {
  e <- numeric(length(lp)); e[k] <- h
  (gfull(lp + e) - gfull(lp - e)) / (2 * h)
})
cat("FD H_rf max:", max(abs(Hfd[r, ])), "\n")
Qg <- as.matrix(sdr$joint_precision_white)
cat("Qg_rf max:", max(abs(Qg[r, f])), " rel diff to FD:",
    max(abs(Qg[r, f] - Hfd[r, ])) / max(abs(Hfd[r, ])), "\n")
cat("Qg_ff:", format(Qg[f, f], digits = 8), "\n")
cat("H_ff FD:", format(Hfd[f, ], digits = 8), "\n")
Hr <- numDeriv::jacobian(gfull, lp, method = "Richardson")
Hr <- (Hr + t(Hr)) / 2
cat("Richardson H_rf vs FD:", max(abs(Hr[r, f] - Hfd[r, ])) / max(abs(Hr[r, f])), "\n")
cat("Richardson H_rf vs Qg_rf:", max(abs(Hr[r, f] - Qg[r, f])) / max(abs(Hr[r, f])), "\n")
cat("Richardson H_rr vs Qg_rr:", max(abs(Hr[r, r] - Qg[r, r])) / max(abs(Hr[r, r])), "\n")
k <- which.max(abs(Hr[r, f] - Qg[r, f]))
cat("worst entry:", arrayInd(k, c(length(r), length(f))), Hr[r, f][k], Qg[r, f][k], "\n")
nm <- names(lp)
D <- abs(Hr[r, f] - Qg[r, f])
cat("by name of r (max abs diff):", tapply(apply(D, 1, max), nm[r], max), "\n")
cat("by name of r (max abs H):", tapply(apply(abs(Hr[r, f]), 1, max), nm[r], max), "\n")
cat("lp == sdr par? ", max(abs(lp[r] - frmtmb:::beta_whiten_map(sdr$par.random, fit$beta_whiten, inv = TRUE)[seq_along(r)])), "\n")
cat("par.random names:", unique(names(sdr$par.random)), "\n")
cat("sdr par.fixed:", format(sdr$par.fixed, digits = 10), " opt:", format(fit$opt$par, digits = 10), "\n")

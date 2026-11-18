# Review of lane sparsese: where the 3.6e-7 of fitted(re_formula = NA) at
# 3510 joint parameters, group sd 3, comes from (both routes agree with
# each other, sparsese-agree.R). The package's whitened precision Qg is
# compared with the analytic one in the package's own coordinates
# (X T_pkg formed first), and the SE of x = 1..2 is recomputed from each.
#   Rscript sparsese-review-floor.R <groups> <seed> <group sd>
a <- commandArgs(trailingOnly = TRUE)
ng <- as.integer(a[1])
seed <- as.integer(a[2])
gsd <- as.numeric(a[3])
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
options(frmtmb.joint_root_max = 100)
set.seed(seed)
n <- 2L * ng
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(ng, 0, gsd)[d$g] + stats::rnorm(n, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fit <- suppressWarnings(suppressMessages(frm(
  y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g), data = d, REML = TRUE)))
sdr <- frmtmb:::sdr_of(fit)
Qg <- sdr$joint_precision_white
r <- fit$obj$env$random
rn <- rownames(Qg)
jb <- which(rn == "beta")
jbb <- which(rn == "b")
Tp <- fit$beta_whiten$blocks[[1]]$T
Xn <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7, d)
XT <- Xn %*% Tp
s2 <- exp(2 * fit$opt$par[["betad"]])
sd2 <- exp(2 * fit$opt$par[["theta"]])
gi <- as.integer(d$g)
Hbb <- crossprod(XT) / s2
Hbz <- t(rowsum(XT, gi, reorder = TRUE)) / s2
Q <- as.matrix(Qg)
cat("rel |Qg_bb - (XT)'(XT)/s2| / max:",
    format(max(abs(Q[jb, jb] - Hbb)) / max(abs(Hbb)), digits = 3), "\n")
cat("elementwise rel max Qg_bb:",
    format(max(abs(Q[jb, jb] - Hbb) / abs(Hbb)), digits = 3), "\n")
cat("rel |Qg_bz - analytic| / max:",
    format(max(abs(Q[jb, jbb] - Hbz)) / max(abs(Hbz)), digits = 3), "\n")
cat("kappa(XT):", format(kappa(XT, exact = TRUE), digits = 3),
    " kappa(X):", format(kappa(Xn, exact = TRUE), digits = 3), "\n")
# conditional SE of the new rows from Qg_rr and from the analytic H_rr
A <- stats::model.matrix(~ x1 + x2 + x3 + x4 + x5 + x6 + x7,
                         data.frame(x1 = 1, x2 = 1, x3 = 1, x4 = 1, x5 = 1,
                                    x6 = 1, x7 = 1)[rep(1, 11), ])
xs <- seq(1, 2, by = 0.1)
A <- cbind(1, sapply(1:7, function(k) xs^k))
B <- A %*% Tp
Hrr <- Q[r, r]
Ha <- Hrr
ib <- match(jb, r)
ibb <- match(jbb, r)
Ha[ib, ib] <- Hbb
Ha[ib, ibb] <- Hbz
Ha[ibb, ib] <- t(Hbz)
Ha[cbind(ibb, ibb)] <- tabulate(gi, ng) / s2 + 1 / sd2
cat("rel |Qg_zz - analytic| (diag):",
    format(max(abs(diag(Hrr)[ibb] - diag(Ha)[ibb]) / diag(Ha)[ibb]),
           digits = 3), "\n")
Wq <- chol2inv(chol(Hrr))[ib, ib]
Wa <- chol2inv(chol(Ha))[ib, ib]
vq <- rowSums((B %*% Wq) * B)
va <- rowSums((B %*% Wa) * B)
cat("conditional SE, Qg vs analytic (package T): rel",
    format(max(abs(sqrt(vq) / sqrt(va) - 1)), digits = 3), "\n")

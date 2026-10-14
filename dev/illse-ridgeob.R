# Lane illse: the nonlinear ridge of test-se-tier.R ("a nonlinear ridge
# with a random effect loses its SEs") under the lane library of install
# 8: what se_reparam_cov() measured in the coordinates of reparam_basis()
# before a design with condition number below 1e4 kept its own
# coordinates. Run under each BLAS.
.libPaths(c("C:/Users/adf44/source/r/wt-illse-lib", "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("R", R.home(), "\n")
set.seed(1); k <- 10
d <- data.frame(f = factor(rep(seq_len(k), each = 6)), g = factor(rep(1:6, k)))
d$y <- rnorm(k)[d$f] + rnorm(6, 0, 0.5)[d$g] + rnorm(k * 6, 0, 0.5)
fit <- suppressWarnings(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1 + (1 | g), nl = TRUE), data = d))
cat("lost:", names(ns$sdr_of(fit)$se_lost), "| root:", !is.null(ns$sdr_of(fit)$cov_fixed_root), "\n")
p <- fit$opt$par
M <- ns$reparam_basis(fit$frame, names(p))
h <- ns$reparam_hessian_refined(fit$obj, p, M)
Dt <- sqrt(diag(h$H)); et <- eigen(h$H / outer(Dt, Dt), TRUE)
nk <- sqrt(colSums(((h$E / outer(Dt, Dt)) %*% et$vectors)^2))
cat("smallest / largest:", format(min(et$values) / max(et$values), digits = 3),
    "| its noise:", format(nk[length(nk)], digits = 3), "\n")

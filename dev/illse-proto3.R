# Lane illse, prototype: gaussian raw polynomials with (1 | g). The outer
# Hessian of the Laplace objective by differences of its gradient along
# reparam_basis(), against the exact marginal covariance at the fit's
# own point (dev/illse-common.R lmm_exact_cov()). Degrees 5 to 8, seeds
# 1 to 3. Also lme4's vcov() over the exact beta block, where lme4 fits.
#   ILLSE_ARM=base Rscript dev/illse-proto3.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
ns <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-illse/R/se-reparam.R",
           envir = ns)
for (deg in 5:8) for (seed in 1:3) {
  ds <- des_poly(deg, re = TRUE)(seed)
  f <- suppressWarnings(frm(ds$formula, data = ds$data))
  p <- f$opt$par
  X <- as.matrix(f$frame$linpreds[[1]]$X)
  Ve <- lmm_exact_cov(X, ds$data$y, ds$data$g, p)
  M <- ns$reparam_basis(f$frame, names(p))
  h <- ns$reparam_hessian(f$obj, p, M)
  D <- sqrt(diag(h$H))
  ev <- eigen(h$H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  V <- M %*% solve(h$H) %*% t(M)
  ee <- mp_cov_err(V, Ve)
  Vo <- suppressWarnings(vcov(f, full = TRUE))
  eo <- if (all(is.finite(Vo))) mp_cov_err(Vo, Ve) else c(cor = NaN, rel = NaN, diag = NaN)
  cat(sprintf("deg %d (1|g) s%d | kappa(M'HM) %.2e | reparam cor %.1e rel %.1e diag %.1e | 0.70.0 vcov(full) cor %.1e diag %.1e\n",
              deg, seed, max(ev) / min(ev), ee[["cor"]], ee[["rel"]],
              ee[["diag"]], eo[["cor"]], eo[["diag"]]))
}

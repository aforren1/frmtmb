# Lane illse, prototype before any change to the package: the outer
# Hessian by central differences of the gradient along the columns of
# reparam_basis() (R/se-reparam.R, sourced into the 0.70.0 namespace),
# inverted there and mapped back, against the exact covariance at the
# fit's own point (256-bit, dev/illse-common.R gauss_exact_cov()).
# Gaussian raw polynomials without random effects, degrees 3 to 8,
# seeds 1 to 3. Reports kappa of M'HM at unit diagonal and the largest
# relative error over all entries of the covariance and its diagonal,
# for the new route and for solve() of the natural AD Hessian.
#   ILLSE_ARM=base Rscript dev/illse-proto.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
ns <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-illse/R/se-reparam.R",
           envir = ns)
for (deg in 3:8) for (seed in 1:3) {
  ds <- des_poly(deg)(seed)
  f <- suppressWarnings(frm(ds$formula, data = ds$data))
  p <- f$opt$par
  X <- f$frame$linpreds[[1]]$X
  Ve <- gauss_exact_cov(as.matrix(X), ds$data$y, p[names(p) == "beta"],
                        p[names(p) == "betad"])
  M <- ns$reparam_basis(f$frame, names(p))
  h <- ns$reparam_hessian(f$obj, p, M)
  D <- sqrt(diag(h$H))
  ev <- eigen(h$H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  C <- solve(h$H)
  V <- M %*% C %*% t(M)
  Hx <- f$obj$he(p)
  Vn <- tryCatch(solve(Hx), error = function(e) matrix(NaN, nrow(Hx), nrow(Hx)))
  en <- if (all(is.finite(Vn))) mp_cov_err(Vn, Ve) else c(cor = NaN, rel = NaN, diag = NaN)
  ee <- mp_cov_err(V, Ve)
  cat(sprintf("deg %d s%d | kappa(M'HM) %.2e | reparam cor %.1e rel %.1e diag %.1e | solve(he) cor %.1e rel %.1e diag %.1e\n",
              deg, seed, max(ev) / min(ev), ee[["cor"]], ee[["rel"]], ee[["diag"]],
              en[["cor"]], en[["rel"]], en[["diag"]]))
}

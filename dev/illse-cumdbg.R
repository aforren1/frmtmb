# Lane illse: which route the covariance of a cumulative raw polynomial
# takes, and the conditioning in each pass (degree from DEG, seed 1).
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
ns <- asNamespace("frmtmb")
ds <- des_poly(as.integer(Sys.getenv("DEG", "7")), FALSE, "cumulative")(1)
f <- suppressWarnings(frm(ds$formula, data = ds$data, family = ds$family))
s <- ns$sdr_of(f)
cat("root:", !is.null(s$cov_fixed_root), " design kappa:",
    ns$se_design_kappa(f), "\n")
p <- f$opt$par
Hx <- f$obj$he(p); D <- sqrt(diag(Hx))
ev <- eigen(Hx / outer(D, D), TRUE, only.values = TRUE)$values
cat("natural kappa:", max(ev) / min(ev), "\n")
M <- ns$reparam_basis(f$frame, names(p))
h1 <- ns$reparam_hessian(f$obj, p, M)
D1 <- sqrt(diag(h1$H)); e1 <- eigen(h1$H / outer(D1, D1), TRUE)$values
cat("first pass kappa:", max(e1) / min(e1), "\n")
h2 <- ns$reparam_hessian_refined(f$obj, p, M)
D2 <- sqrt(diag(h2$H)); e2 <- eigen(h2$H / outer(D2, D2), TRUE)$values
cat("refined kappa:", max(e2) / min(e2), "\n")

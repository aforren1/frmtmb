# Reviewer of lane seflags, check of punch 2b: a bernoulli fit of the
# degree-8 raw polynomial (seed 8, median split) whose direction stays
# "flat": the null basis, the flat subspace's tolerance, the rows'
# cosines with it, and the rows' SEs.
#   Rscript dev/seflags-rev4-flatglm.R <lane|base>

source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev4-common.R")
set.seed(8)
d <- data.frame(x = stats::runif(200, 1, 2))
d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
d$yb <- as.integer(d$y > stats::median(d$y))
rhs <- paste0("x", 1:8, collapse = " + ")
f <- suppressWarnings(frm(as.formula(paste("yb ~", rhs)), data = d, family = bernoulli()))
sdr <- suppressWarnings(ns$sdr_of(f))
cat("se_lost:", names(sdr$se_lost), "\n")
cat("se_null:", if (is.null(sdr$se_null)) "NULL" else dim(sdr$se_null), "\n")
fl <- sdr$se_null_flat
cat("flat cols:", fl$cols, "tol", fl$tol, "\n")
cat("se_ill:", if (is.null(sdr$se_ill)) "NULL" else length(sdr$se_ill$lam), "lam", sdr$se_ill$lam, "\n")
print(signif(sdr$se_null, 3))
H <- f$obj$he(f$opt$par); D <- sqrt(diag(H)); ev <- eigen(H / outer(D, D), symmetric = TRUE)$values
cat("unit-diag eigenvalues / max:", signif(ev / max(ev), 3), "\n")
nd <- data.frame(x = c(1.5, 2, 2.5)); for (k in 1:8) nd[[paste0("x", k)]] <- nd$x^k
X0 <- model.matrix(as.formula(paste("~", rhs)), nd)
N <- sdr$se_null
if (!is.null(N) && NCOL(N)) {
  s <- fl$dscale; s[!is.finite(s) | !(s > 0)] <- 1
  Nq <- N * s; Gq <- X0 / rep(s, each = 3)
  cat("cosines of the rows with the null columns (unit-diagonal):\n")
  print(signif(abs(Gq %*% Nq) / outer(sqrt(rowSums(Gq^2)), sqrt(colSums(Nq^2))), 3))
}
lb <- frm_lp_basis(f, newdata = nd)
cat("se_nonest:", lb$se_nonest, "| A V A' SE:", signif(sqrt(rowSums((as.matrix(lb$A) %*% lb$V) * as.matrix(lb$A)) + lb$extra_var), 4), "\n")

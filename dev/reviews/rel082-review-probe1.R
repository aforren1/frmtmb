# Release review of 0.82.0, probe 1: brms's pars rule, Matrix's Cholesky
# on a matrix that is not positive definite, TMB's spHess layout.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")), "\n")
cat("brms", as.character(packageVersion("brms")), "Matrix",
    as.character(packageVersion("Matrix")), "\n")
print(brms:::extract_pars)

# Cholesky of a sparse matrix that is not PD: warning or error?
library(Matrix)
A <- Matrix(c(2, 1, 0, 1, -1, 0, 0, 0, 3), 3, sparse = TRUE)
A <- forceSymmetric(A)
r <- tryCatch(withCallingHandlers({
  Cholesky(A, LDL = FALSE, perm = TRUE); "no condition"
}, warning = function(w) {
  cat("WARNING:", conditionMessage(w), "\n")
}), error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
# near-singular, PD: 1e-14 on the diagonal after the shift
B <- Diagonal(3, 1e-14) + 0
r2 <- tryCatch(withCallingHandlers({
  Cholesky(forceSymmetric(as(B, "CsparseMatrix")), LDL = FALSE, perm = TRUE)
  "no condition"
}, warning = function(w) cat("WARNING2:", conditionMessage(w), "\n")),
error = function(e) paste("ERROR2:", conditionMessage(e)))
print(r2)
# a large one, supernodal
set.seed(1)
n <- 3000
M <- rsparsematrix(n, n, 0.002)
S <- forceSymmetric(crossprod(M) + Diagonal(n, 1e-3))
S2 <- S - Diagonal(n, 10)
r3 <- tryCatch(withCallingHandlers({
  Cholesky(S2, LDL = FALSE, perm = TRUE); "no condition"
}, warning = function(w) cat("WARNING3:", conditionMessage(w), "\n")),
error = function(e) paste("ERROR3:", conditionMessage(e)))
print(r3)
# with NaN
S3 <- S
S3[1, 1] <- NaN
r4 <- tryCatch(withCallingHandlers({
  Cholesky(S3, LDL = FALSE, perm = TRUE); "no condition"
}, warning = function(w) cat("WARNING4:", conditionMessage(w), "\n")),
error = function(e) paste("ERROR4:", conditionMessage(e)))
print(r4)

# TMB's spHess(random = TRUE): class and triangle
library(frmtmb)
set.seed(2)
d <- data.frame(g = factor(rep(1:30, each = 5)), x = rnorm(150))
d$y <- rbinom(150, 1, plogis(0.3 * d$x + rnorm(30)[d$g]))
f <- frm(y ~ x + (1 + x | g), data = d, family = bernoulli())
env <- f$obj$env
f$obj$fn(f$opt$par)
H <- env$spHess(env$last.par, random = TRUE)
cat("spHess class:", class(H), "\n")
if (methods::is(H, "symmetricMatrix")) cat("uplo:", H@uplo, "\n")
Hg <- as(H, "CsparseMatrix")
A1 <- forceSymmetric(Hg)
full <- as.matrix(A1)
cat("off-diagonal nonzeros kept by forceSymmetric:", sum(full[upper.tri(full)] != 0),
    " lower in H:", sum(as.matrix(Hg)[lower.tri(full)] != 0),
    " upper in H:", sum(as.matrix(Hg)[upper.tri(full)] != 0), "\n")
cat("min eigen (forceSymmetric):", min(eigen(full, TRUE, only.values = TRUE)$values), "\n")
# against the dense full Hessian
Hd <- as.matrix(Hg)
Hd <- Hd + t(Hd) - diag(diag(Hd))
cat("min eigen (symmetrized from lower):",
    min(eigen(Hd, TRUE, only.values = TRUE)$values), "\n")
cat("gp_inner_pd at 1e-2:", frmtmb:::gp_inner_pd(f$obj, f$opt$par, 1e-2), "\n")

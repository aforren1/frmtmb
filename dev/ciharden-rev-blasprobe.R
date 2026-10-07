# Reviewer probe of one R build: which BLAS and LAPACK it runs, how many
# OpenBLAS threads, and a fingerprint of BLAS-only and LAPACK-only
# results, so builds can be compared bit for bit.
#   Rscript dev/ciharden-rev-blasprobe.R <probe dll>
a <- commandArgs(trailingOnly = TRUE)
dyn.load(a[1])
cat("R.home:", R.home(), "\n")
cat("La_library:", La_library(), "\n")
cat("La_version:", La_version(), "\n")
esv <- extSoftVersion()
cat("extSoftVersion BLAS:", esv[["BLAS"]], "\n")
cat("sessionInfo BLAS/LAPACK:", sessionInfo()$BLAS, "|",
    sessionInfo()$LAPACK, "\n")
cat("env OPENBLAS_NUM_THREADS:", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
info <- .Call("rev_ob_info")
cat("openblas dll:", info[1], "\n")
cat("openblas_get_num_threads:", info[2], "\n")
cat("openblas_get_config:", info[3], "\n")
hx <- function(x) paste(sprintf("%a", sum(x * seq_along(x))), collapse = "")
set.seed(1)
A <- matrix(rnorm(400 * 400), 400)
S <- crossprod(A) + diag(400)
B <- matrix(rnorm(400 * 3), 400)
# BLAS level 3 only
cat("fp dgemm:", hx(A %*% A), "\n")
# LAPACK: Cholesky, LU solve, QR (LAPACK = TRUE), symmetric eigen, SVD
cat("fp dpotrf:", hx(chol(S)), "\n")
cat("fp dgesv:", hx(solve(A, B)), "\n")
cat("fp dgeqp3:", hx(qr.R(qr(A, LAPACK = TRUE))), "\n")
cat("fp dsyevr:", hx(eigen(S, symmetric = TRUE, only.values = TRUE)$values),
    "\n")
cat("fp dgesdd:", hx(svd(A, nu = 0, nv = 0)$d), "\n")
# timing, which a threaded OpenBLAS changes most for big LAPACK calls
S2 <- crossprod(matrix(rnorm(3000 * 3000), 3000)) + diag(3000)
cat("t dpotrf3000:", system.time(chol(S2))[["elapsed"]], "\n")
G <- matrix(rnorm(3000 * 3000), 3000)
cat("t dgemm3000:", system.time(G %*% G)[["elapsed"]], "\n")

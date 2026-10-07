# Is this R's BLAS (and LAPACK) OpenBLAS? A 2000 x 2000 product takes
# about 4 s with the reference BLAS and 0.05 s with OpenBLAS; a Cholesky
# of a 3000 x 3000 matrix through LAPACK's dpotrf about 1.5 s with the
# reference LAPACK and well under that with OpenBLAS's own dpotrf.
# OPENBLAS_NUM_THREADS is printed because it decides the rounding of
# threaded kernels.
cat("R:", R.home(), "\n")
cat("OPENBLAS_NUM_THREADS:", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
print(La_version())
set.seed(1)
m <- matrix(rnorm(4e6), 2000)
cat("dgemm 2000 (s):", system.time(m %*% m)[["elapsed"]], "\n")
s <- crossprod(matrix(rnorm(9e6), 3000)) + diag(3000)
cat("dpotrf 3000 (s):", system.time(chol(s))[["elapsed"]], "\n")
cat("dgesv 2000 (s):", system.time(solve(m))[["elapsed"]], "\n")

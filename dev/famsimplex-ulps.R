# The density tests' error measure after fix round 1: |ours - reference|
# over the reference's sum of absolute terms, in double epsilons. On the
# tests' own draws (seeds 4, 5, 6 of test-simplex-families.R) and on 400
# wider random cases (seed 20261014, K 2 to 6, phi to exp(8)).
# Output: dev/famsimplex-ulps.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
ULP <- .Machine$double.eps
sm <- function(E) { E <- E - apply(E, 1, max); P <- exp(E); P / rowSums(P) }
rd <- function(A) { G <- matrix(rgamma(length(A), A), nrow(A)); G / rowSums(G) }
ref_d <- function(Y, A) {
  t1 <- lgamma(rowSums(A)); t2 <- lgamma(A); t3 <- (A - 1) * log(Y)
  list(v = t1 - rowSums(t2) + rowSums(t3),
       s = abs(t1) + rowSums(abs(t2)) + rowSums(abs(t3)))
}
ref_l <- function(Y, M, S, C, ref) {
  X <- log(Y[, -ref, drop = FALSE]) - log(Y[, ref]); D <- ncol(X)
  p <- vapply(seq_len(nrow(Y)), function(i) {
    Sig <- diag(S[i, ], D) %*% C %*% diag(S[i, ], D); r <- X[i, ] - M[i, ]
    c(D * log(2 * pi), as.numeric(determinant(Sig)$modulus),
      sum(r * solve(Sig, r)), sum(log(Y[i, ])), sum(abs(log(Y[i, ]))))
  }, numeric(5))
  list(v = -0.5 * colSums(p[1:3, , drop = FALSE]) - p[4, ],
       s = 0.5 * colSums(abs(p[1:3, , drop = FALSE])) + p[5, ])
}
run_d <- function(K, ref, n, a, phimax, sdE) {
  cats <- letters[seq_len(K)]
  fam <- ns$fam_dirichlet_impl(cats, ref); sx <- fam$simplex
  Y <- rd(matrix(a, n, K))
  E <- matrix(rnorm(n * length(sx$pred), 0, sdE), n)
  phi <- exp(runif(n, -1, phimax)); M <- matrix(0, n, K); M[, sx$pred] <- E
  dp <- c(setNames(lapply(seq_along(sx$pred), function(j) E[, j]), sx$mu),
          list(phi = phi))
  o <- fam$lpdf(Y, dp, list()); r <- ref_d(Y, sm(M) * phi)
  c(scaled = max(abs(o - r$v) / r$s) / ULP,
    old = max(abs(o - r$v) / pmax(1, abs(r$v))) / ULP)
}
run_l <- function(K, ref, n) {
  cats <- letters[seq_len(K)]; D <- K - 1L
  fam <- ns$fam_logistic_normal_impl(cats, ref); sx <- fam$simplex
  Y <- rd(matrix(2, n, K)); M <- matrix(rnorm(n * D), n)
  S <- matrix(exp(runif(n * D, -1, 1)), n)
  C <- if (D > 1) cov2cor(crossprod(matrix(rnorm(D * (D + 3)), D + 3))) else
    matrix(1, 1, 1)
  dp <- c(setNames(lapply(seq_len(D), function(j) M[, j]), sx$mu),
          setNames(lapply(seq_len(D), function(j) S[, j]), sx$sigma))
  o <- if (D > 1) fam$lpdf(Y, dp, list(), list(lncor_raw = ns$us_theta_cor(C)))
       else fam$lpdf(Y, dp, list())
  r <- ref_l(Y, M, S, C, ref)
  c(scaled = max(abs(o - r$v) / r$s) / ULP,
    old = max(abs(o - r$v) / pmax(1, abs(r$v))) / ULP)
}
cat("\n== the tests' own draws\n")
set.seed(4)
td <- do.call(rbind, lapply(2:4, function(K) do.call(rbind, lapply(
  c(1L, K, NA_integer_), function(ref) run_d(K, ref, 50, 1.5, 4, 1.5)))))
cat("dirichlet  worst scaled", round(max(td[, "scaled"]), 2),
    " worst old", round(max(td[, "old"]), 1), "\n")
set.seed(6)
tl <- do.call(rbind, lapply(2:4, function(K) do.call(rbind, lapply(
  unique(c(1L, K)), function(ref) run_l(K, ref, 40)))))
cat("logistic_normal worst scaled", round(max(tl[, "scaled"]), 2),
    " worst old", round(max(tl[, "old"]), 1), "\n")
cat("\n== 400 wider random cases\n")
set.seed(20261014)
wd <- t(replicate(400, {
  K <- sample(2:6, 1)
  run_d(K, sample(c(seq_len(K), NA_integer_), 1), 20, runif(1, 0.2, 5),
        8, 3)
}))
cat("dirichlet  worst scaled", round(max(wd[, "scaled"]), 2),
    " worst old", round(max(wd[, "old"]), 1), "\n")
wl <- t(replicate(400, { K <- sample(2:6, 1); run_l(K, sample(K, 1), 20) }))
cat("logistic_normal worst scaled", round(max(wl[, "scaled"]), 2),
    " worst old", round(max(wl[, "old"]), 1), "\n")
cat("test bound: 64\n")

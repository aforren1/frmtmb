# Review of lane arspat2: does determinant(G, logarithm = TRUE) on an RTMB
# advector matrix (autocor_seam_t_margin(), R/autocor-seam.R) survive a
# determinant outside double range? G = s * I_m and G = s * Q (Q a fixed
# SPD matrix), log|G| = m log s + log|Q|, taped through MakeADFun.
source("dev/reviews/arspat2-review-common.R")
cat("RTMB", format(packageVersion("RTMB")), "\n")
set.seed(1)
for (m in c(50, 200, 300)) {
  Q <- crossprod(matrix(rnorm(m * m), m)) / m + diag(m)
  ldQ <- as.numeric(determinant(Q)$modulus)
  for (s in c(1, 0.05, 20)) {
    f <- function(p) {
      G <- RTMB::matrix(exp(p$ls) * as.vector(Q), m, m)
      determinant(G, logarithm = TRUE)[["modulus"]]
    }
    obj <- RTMB::MakeADFun(f, list(ls = log(s)), silent = TRUE)
    v <- obj$fn(log(s))
    g <- obj$gr(log(s))
    cat(sprintf("m %3d s %5.2f: taped log|G| %s want %.4f, gradient %s want %d\n",
                m, s, format(v), m * log(s) + ldQ, format(g), m))
  }
}

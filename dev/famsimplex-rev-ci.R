# Reviewer: confint rows of family extras vs rescor; hypothesis profile on
# lncor; the margins of the lane's density tests. Output:
# dev/famsimplex-rev-ci.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
show <- function(lab, expr) {
  cat("\n---", lab, "\n")
  r <- tryCatch(suppressWarnings(expr),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r) && length(r) == 1) cat(r, "\n") else print(r)
}
set.seed(2)
d <- data.frame(x = rnorm(200))
e <- MASS::mvrnorm(200, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2))
d$y1 <- 1 + d$x + e[, 1]
d$y2 <- 2 - d$x + e[, 2]
fr <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE), data = d)
show("rescor confint rows", tail(confint(fr), 3))
show("rescor confint parm rescor__y1__y2",
     confint(fr, parm = "rescor__y1__y2"))
fl <- frm(bf(Y ~ x), family = logistic_normal(),
          data = fs_data_logistic_normal(seed = 13))
show("ln hypothesis profile", hypothesis(fl, "lncor__b__c = 0",
                                         class = NULL, method = "profile"))
# margins of the lane's density tests, same seeds and draws as the tests
ULP <- .Machine$double.eps
sx_softmax <- fs_softmax
sx_rdirichlet <- fs_rdirichlet
ref_dirichlet <- function(Y, A) {
  lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A - 1) * log(Y))
}
set.seed(4)
m <- 0
for (K in 2:4) {
  cats <- letters[seq_len(K)]
  for (ref in c(1L, K, NA_integer_)) {
    fam <- frmtmb:::fam_dirichlet_impl(cats, ref)
    sx <- fam$simplex
    n <- 50
    Y <- sx_rdirichlet(matrix(1.5, n, K))
    E <- matrix(stats::rnorm(n * length(sx$pred), 0, 1.5), n)
    phi <- exp(stats::runif(n, -1, 4))
    M <- matrix(0, n, K)
    M[, sx$pred] <- E
    dp <- c(stats::setNames(lapply(seq_along(sx$pred),
                                   function(j) E[, j]), sx$mu),
            list(phi = phi))
    ours <- fam$lpdf(Y, dp, list())
    ref_v <- ref_dirichlet(Y, sx_softmax(M) * phi)
    m <- max(m, max(abs(ours - ref_v) / pmax(1, abs(ref_v))))
  }
}
cat("test-simplex-families dirichlet density: worst", signif(m / ULP, 4),
    "ulps against a bound of 1000\n")
set.seed(5)
fam <- frmtmb:::fam_dirichlet2_impl(c("a", "b", "c"))
Y <- sx_rdirichlet(matrix(2, 60, 3))
A <- matrix(exp(stats::runif(180, -2, 4)), 60)
dp <- stats::setNames(lapply(1:3, function(k) A[, k]), fam$simplex$mu)
m2 <- max(abs(fam$lpdf(Y, dp, list()) - ref_dirichlet(Y, A)) /
            pmax(1, abs(ref_dirichlet(Y, A))))
cat("dirichlet2 density: worst", signif(m2 / ULP, 4), "ulps\n")
# the same check with a condition-aware scale: the sum of |terms|
set.seed(20261013)
w1 <- w2 <- 0
for (r in 1:200) {
  K <- sample(2:6, 1)
  n <- 30
  Y <- fs_rdirichlet(matrix(exp(runif(n * K, -1, 2)), n))
  fam <- frmtmb:::fam_dirichlet_impl(paste0("k", 1:K), 1L)
  E <- matrix(rnorm(n * (K - 1), 0, 3), n)
  phi <- exp(runif(n, -2, 6))
  dp <- c(setNames(lapply(1:(K - 1), function(j) E[, j]), fam$simplex$mu),
          list(phi = phi))
  A <- fs_softmax(cbind(0, E)) * phi
  refv <- ref_dirichlet(Y, A)
  scale <- abs(lgamma(rowSums(A))) + rowSums(abs(lgamma(A))) +
    rowSums(abs((A - 1) * log(Y)))
  ours <- fam$lpdf(Y, dp, list())
  w1 <- max(w1, max(abs(ours - refv) / pmax(1, abs(refv))))
  w2 <- max(w2, max(abs(ours - refv) / scale))
}
cat("200 random dirichlet cases: worst relative to max(1,|ll|)",
    signif(w1 / ULP, 4), "ulps; relative to the sum of |terms|",
    signif(w2 / ULP, 4), "ulps\n")

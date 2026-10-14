# Reviewer: densities against independent references (mvtnorm, lgamma),
# at random points and at fitted estimates through the whole objective.
# Output: dev/famsimplex-rev-dens.log
.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
set.seed(20261013)
worst <- c(dir = 0, ln = 0, dir2 = 0)
for (rep in 1:40) {
  K <- sample(2:6, 1)
  cats <- paste0("k", seq_len(K))
  n <- 30
  Y <- fs_rdirichlet(matrix(exp(runif(n * K, -1, 2)), n))
  # dirichlet, a random reference or NA
  ref <- sample(c(seq_len(K), NA), 1)
  fam <- frmtmb:::fam_dirichlet_impl(cats, ref)
  pr <- fam$simplex$pred
  E <- matrix(rnorm(n * length(pr), 0, 3), n)
  phi <- exp(runif(n, -2, 6))
  dp <- c(setNames(lapply(seq_along(pr), function(j) E[, j]),
                   fam$simplex$mu), list(phi = phi))
  M <- matrix(0, n, K)
  M[, pr] <- E
  A <- fs_softmax(M) * phi
  refv <- lgamma(rowSums(A)) - rowSums(lgamma(A)) + rowSums((A - 1) * log(Y))
  ours <- fam$lpdf(Y, dp, list())
  worst["dir"] <- max(worst["dir"], max(abs(ours - refv) / pmax(1, abs(refv))))
  # dirichlet2
  fam2 <- frmtmb:::fam_dirichlet2_impl(cats)
  A2 <- matrix(exp(runif(n * K, -3, 5)), n)
  dp2 <- setNames(lapply(seq_len(K), function(k) A2[, k]), fam2$simplex$mu)
  refv2 <- lgamma(rowSums(A2)) - rowSums(lgamma(A2)) +
    rowSums((A2 - 1) * log(Y))
  ours2 <- fam2$lpdf(Y, dp2, list())
  worst["dir2"] <- max(worst["dir2"],
                       max(abs(ours2 - refv2) / pmax(1, abs(refv2))))
  # logistic_normal against mvtnorm::dmvnorm on the alr scale + Jacobian
  rl <- sample(seq_len(K), 1)
  faml <- frmtmb:::fam_logistic_normal_impl(cats, rl)
  D <- K - 1
  pl <- faml$simplex$pred
  Mu <- matrix(rnorm(n * D), n)
  S <- matrix(exp(runif(n * D, -1.5, 1.5)), n)
  C <- if (D > 1) cov2cor(crossprod(matrix(rnorm(D * (D + 2)), D + 2))) else
    matrix(1, 1, 1)
  dpl <- c(setNames(lapply(seq_len(D), function(j) Mu[, j]), faml$simplex$mu),
           setNames(lapply(seq_len(D), function(j) S[, j]),
                    faml$simplex$sigma))
  X <- log(Y[, pl, drop = FALSE]) - log(Y[, rl])
  refl <- vapply(seq_len(n), function(i) {
    Sig <- diag(S[i, ], D) %*% C %*% diag(S[i, ], D)
    mvtnorm::dmvnorm(X[i, ], Mu[i, ], Sig, log = TRUE) - sum(log(Y[i, ]))
  }, 0)
  oursl <- if (D > 1) {
    faml$lpdf(Y, dpl, list(), list(lncor_raw = frmtmb:::us_theta_cor(C)))
  } else faml$lpdf(Y, dpl, list())
  worst["ln"] <- max(worst["ln"], max(abs(oursl - refl) / pmax(1, abs(refl))))
}
cat("worst relative error over 40 random cases (K 2..6, n 30):\n")
print(signif(worst, 3))
cat("in ulps:", signif(worst / .Machine$double.eps, 3), "\n")
# whole objective at the estimates: K = 4 ln, refcat c, sigma<k> ~ z
set.seed(41)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
Mu <- cbind(0.2 + 0.4 * d$x, -0.1 * d$x, 0.5)
C <- matrix(c(1, 0.3, -0.2, 0.3, 1, 0.4, -0.2, 0.4, 1), 3)
Sd <- cbind(exp(-0.7 + 0.3 * d$z), 0.7, 0.9)
Y <- t(vapply(seq_len(n), function(i) {
  fs_rlogistic_normal(Mu[i, , drop = FALSE],
                      diag(Sd[i, ]) %*% C %*% diag(Sd[i, ]), ref = 3)[1, ]
}, numeric(4)))
colnames(Y) <- c("a", "b", "c", "d")
d$Y <- Y
f <- frm(bf(Y ~ x, sigmaa ~ z), family = logistic_normal(refcat = "c"),
         data = d)
b <- fixef(f)[, "Estimate"]
sp <- summary(f)$spec_pars
Mh <- cbind(b[["mua_Intercept"]] + b[["mua_x"]] * d$x,
            b[["mub_Intercept"]] + b[["mub_x"]] * d$x,
            b[["mud_Intercept"]] + b[["mud_x"]] * d$x)
Sh <- cbind(exp(b[["sigmaa_Intercept"]] + b[["sigmaa_z"]] * d$z),
            sp["sigmab", "Estimate"], sp["sigmad", "Estimate"])
Ch <- diag(3)
Ch[1, 2] <- Ch[2, 1] <- sp["lncor(a,b)", "Estimate"]
Ch[1, 3] <- Ch[3, 1] <- sp["lncor(a,d)", "Estimate"]
Ch[2, 3] <- Ch[3, 2] <- sp["lncor(b,d)", "Estimate"]
X <- log(Y[, c(1, 2, 4)]) - log(Y[, 3])
ll <- sum(vapply(seq_len(n), function(i) {
  mvtnorm::dmvnorm(X[i, ], Mh[i, ], diag(Sh[i, ]) %*% Ch %*% diag(Sh[i, ]),
                   log = TRUE) - sum(log(Y[i, ]))
}, 0))
cat(sprintf("ln K=4 refcat c: logLik %.10f mvtnorm %.10f rel %.3e\n",
            logLik(f), ll, abs(logLik(f) - ll) / abs(ll)))
print(round(sp, 4))
cat("truth: sigmab 0.7 sigmad 0.9; cor(a,b) 0.3 (a,d) -0.2 (b,d) 0.4\n")

# REVIEW script 01, claim 2: the Student-t branch of
# rescor_row_loglik() across nu, against THIRD-PARTY references.
#
# Two independent references, because neither covers the whole range:
#   * mvtnorm::dmvt(log = TRUE) is the multivariate t with scale matrix
#     `sigma`, which is what the branch computes. It is built on
#     lgamma() itself, so it is a reference only where lgamma()'s
#     cancellation is harmless: nu up to about 1e6.
#   * mvtnorm::dmvnorm is the nu -> Inf limit, and the t approaches it
#     as O(1/nu), so at nu = 1e8 and above the limit is the reference
#     and the residual must FALL like 1/nu rather than grow.
# Both are outside frmtmb, so agreement is a MEASUREMENT.
#
#   Rscript dev/arcovsample-rev-01-tnu.R

LIB <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
stopifnot("arma_cond_resp" %in% getNamespaceExports("frmtmb"))
cat("build: worker lane lib (arma_cond_resp present)\n")
cat("mvtnorm ", format(packageVersion("mvtnorm")), "\n\n")

lsd <- frmtmb:::lgamma_shift_diff
mvtstd <- frmtmb:::mvt_std_loglik

# ---- Part A: the arithmetic, away from any fit ------------------------
# K = 3 so that K/2 is not an integer multiple of 1/2 that could hide a
# shift bug, n = 6 rows, a correlation matrix that is not the identity.
set.seed(99L)
K <- 3L
n <- 6L
A <- matrix(rnorm(K * K), K)
C <- cov2cor(A %*% t(A) + diag(K))
sg <- c(0.7, 1.3, 2.1)
mu <- matrix(rnorm(n * K), n, K)
y <- mu + matrix(rnorm(n * K), n, K) %*% chol(C) * rep(sg, each = n)
Z <- (y - mu) / rep(sg, each = n)
lsig <- rep(sum(log(sg)), n)
q <- rowSums((Z %*% solve(C)) * Z)
ldet <- as.numeric(determinant(C, logarithm = TRUE)$modulus)
Sigma <- diag(sg) %*% C %*% diag(sg)

new_form <- function(nu) {
  lsd(nu / 2, K / 2) - K / 2 * (log(nu) + log(pi)) - ldet / 2 -
    (nu + K) / 2 * log1p(q / nu) - lsig
}
old_form <- function(nu) {
  lgamma((nu + K) / 2) - lgamma(nu / 2) - K / 2 * log(nu * pi) -
    ldet / 2 - (nu + K) / 2 * log1p(q / nu) - lsig
}
obj_form <- function(nu) {
  vapply(seq_len(n), function(i) {
    as.numeric(mvtstd(matrix(Z[i, ], 1L), C, nu))
  }, numeric(1)) - lsig
}
ref_t <- function(nu) {
  vapply(seq_len(n), function(i) {
    mvtnorm::dmvt(y[i, ], delta = mu[i, ], sigma = Sigma, df = nu,
                  log = TRUE)
  }, numeric(1))
}
ref_norm <- vapply(seq_len(n), function(i) {
  mvtnorm::dmvnorm(y[i, ], mean = mu[i, ], sigma = Sigma, log = TRUE)
}, numeric(1))

nus <- c(1.5, 3, 10, 100, 1e4, 1e6, 1e8, 1e12, 1e20, 1e306)
cat("A. K = 3, n = 6, seed 99. Column meanings:\n")
cat("   dmvt   : max |new - mvtnorm::dmvt|      (valid to ~1e6)\n")
cat("   norm   : max |new - mvtnorm::dmvnorm|   (the nu -> Inf limit)\n")
cat("   objf   : max |new - mvt_std_loglik|     (frmtmb's own tape form)\n")
cat("   oldref : max |OLD - the better reference|\n\n")
cat(sprintf("%9s %14s %14s %14s %14s\n", "nu", "dmvt", "norm", "objf",
            "oldref"))
for (nu in nus) {
  nw <- new_form(nu)
  od <- old_form(nu)
  ob <- obj_form(nu)
  rt <- suppressWarnings(tryCatch(ref_t(nu), error = function(e) rep(NA, n)))
  d_t <- max(abs(nw - rt))
  d_n <- max(abs(nw - ref_norm))
  d_o <- max(abs(nw - ob))
  best <- if (nu <= 1e6) rt else ref_norm
  d_old <- max(abs(od - best))
  cat(sprintf("%9.1e %14.6g %14.6g %14.6g %14.6g\n", nu, d_t, d_n, d_o,
              d_old))
}

cat("\nA2. the SMALL-nu end at full precision, new form vs mvtnorm::dmvt\n")
cat(sprintf("%9s %16s %16s %14s\n", "nu", "new sum", "dmvt sum", "rel"))
for (nu in c(1.05, 1.5, 2, 2.5, 3, 4, 7, 10, 30, 50)) {
  nw <- sum(new_form(nu)); rt <- sum(ref_t(nu))
  cat(sprintf("%9.4g %16.10f %16.10f %14.4g\n", nu, nw, rt,
              abs(nw - rt) / abs(rt)))
}

cat("\nA3. does the residual against the gaussian limit fall like 1/nu?\n")
for (nu in c(1e4, 1e5, 1e6, 1e7, 1e8, 1e10)) {
  cat(sprintf("  nu %8.0e  max|new - dmvnorm| %12.6g  times nu %10.4g\n",
              nu, max(abs(new_form(nu) - ref_norm)),
              max(abs(new_form(nu) - ref_norm)) * nu))
}

# ---- Part B: through the exported function, on a real fit -------------
# Own design and seed, K = 2 with a distributional sigma on one
# response, so the branch is reached with a row-varying sigma.
cat("\nB. rescor_row_loglik() on a fit, nu set through its logm1 link\n")
set.seed(4242L)
N <- 60L
xx <- rnorm(N); zz <- runif(N, -1, 1)
E <- matrix(rnorm(2 * N), N) %*% chol(matrix(c(1, -0.45, -0.45, 1), 2)) /
  sqrt(rchisq(N, 6) / 6)
dd <- data.frame(x = xx, z = zz, y1 = 0.5 + 1.2 * xx + E[, 1],
                 y2 = -0.2 * xx + exp(0.25 * zz) * E[, 2])
fit <- frm(bf(y1 ~ x) + bf(y2 ~ x, sigma ~ z) + set_rescor(TRUE) +
             student(), data = dd)
cat("  logLik at the optimum: ", format(as.numeric(logLik(fit)),
                                        digits = 10), "\n")
bn <- names(fit$frame$par_template$betad)
jnu <- match("nu_(Intercept)", bn)
stopifnot(!is.na(jnu))
Cm <- frmtmb:::us_chol_cor(fit$estimates[["thetar"]], 2L)

cat(sprintf("\n%9s %16s %16s %14s %14s %14s\n", "nu", "sum rrl",
            "-obj$fn", "obj resid", "vs dmvt", "vs dmvnorm"))
for (nu in c(1.5, 3, 10, 100, 1e4, 1e6, 1e8, 1e12, 1e20, 1e306)) {
  p <- fit$opt$par
  p[which(names(p) == "betad")[jnu]] <- log(nu - 1)
  objv <- -as.numeric(fit$obj$fn(p))
  f2 <- fit
  f2$estimates <- fit$obj$env$parList(p)
  f2$cache <- new.env(parent = emptyenv())
  dpv <- frmtmb::eval_dpars(f2)
  rl <- frmtmb::rescor_row_loglik(f2, dpv)
  C2 <- frmtmb:::us_chol_cor(f2$estimates[["thetar"]], 2L)
  m1 <- as.numeric(dpv$y1$mu); m2 <- as.numeric(dpv$y2$mu)
  s1 <- rep_len(as.numeric(dpv$y1$sigma), N)
  s2 <- rep_len(as.numeric(dpv$y2$sigma), N)
  rt <- vapply(seq_len(N), function(i) {
    S <- diag(c(s1[i], s2[i])) %*% C2 %*% diag(c(s1[i], s2[i]))
    suppressWarnings(mvtnorm::dmvt(c(dd$y1[i], dd$y2[i]),
                                   delta = c(m1[i], m2[i]), sigma = S,
                                   df = nu, log = TRUE))
  }, numeric(1))
  rn <- vapply(seq_len(N), function(i) {
    S <- diag(c(s1[i], s2[i])) %*% C2 %*% diag(c(s1[i], s2[i]))
    mvtnorm::dmvnorm(c(dd$y1[i], dd$y2[i]), mean = c(m1[i], m2[i]),
                     sigma = S, log = TRUE)
  }, numeric(1))
  cat(sprintf("%9.1e %16.8f %16.8f %14.6g %14.6g %14.6g\n", nu, sum(rl),
              objv, abs(sum(rl) - objv), max(abs(rl - rt)),
              max(abs(rl - rn))))
}
cat("\nDONE\n")

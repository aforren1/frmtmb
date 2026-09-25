# Lane mv, 2026-09-25: student() with rescor = TRUE.
#
# 1. The taped objective at a SHARED parameter point against
#    mvtnorm::dmvt (an independent implementation of the multivariate t)
#    and against rescor_row_loglik() (plain R, brms's
#    multi_student_t_lpdf term by term).
# 2. The ML fit against a hand-written RTMB objective of the same
#    density, optimized separately.
# 3. The limit nu -> Inf against the gaussian rescor objective.
#
# Run: Rscript dev/mv-validate-student.R  (lane library first on
# .libPaths()). Seeds are fixed below.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-mv")
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))

set.seed(20260925)
n <- 400
K <- 3
x <- rnorm(n)
z <- rnorm(n)
R0 <- matrix(c(1, 0.5, 0.3, 0.5, 1, -0.2, 0.3, -0.2, 1), 3)
W <- rchisq(n, 5) / 5
E <- (matrix(rnorm(n * K), n) %*% chol(R0)) / sqrt(W)
dd <- data.frame(x = x, z = z,
                 y1 = 1 + 0.5 * x + 0.8 * E[, 1],
                 y2 = -1 + 0.3 * x + exp(0.2 * z) * E[, 2],
                 y3 = 0.2 * x + 1.2 * E[, 3])

form <- bf(y1 ~ x) + bf(y2 ~ x, sigma ~ z) + bf(y3 ~ x) +
  set_rescor(TRUE) + student()
fit <- frm(form, data = dd)

## 1. shared parameter point: perturb the optimum so the check is not
## only at a stationary point
obj <- fit$obj
set.seed(7)
p <- fit$opt$par + rnorm(length(fit$opt$par), 0, 0.05)
nll_frm <- obj$fn(p)
pl <- obj$env$parList(p)
beta <- pl$beta
betad <- pl$betad
# layout, read off the template names
bn <- names(fit$frame$par_template$beta)
dn <- names(fit$frame$par_template$betad)
print(bn)
print(dn)
mu1 <- beta[1] + beta[2] * x
mu2 <- beta[3] + beta[4] * x
mu3 <- beta[5] + beta[6] * x
s1 <- exp(betad[match("y1_sigma_(Intercept)", dn)])
s2 <- exp(betad[match("y2_sigma_(Intercept)", dn)] +
            betad[match("y2_sigma_z", dn)] * z)
s3 <- exp(betad[match("y3_sigma_(Intercept)", dn)])
nu <- 1 + exp(betad[match("nu_(Intercept)", dn)])
C <- frmtmb::us_chol_cor(pl$thetar, K)
Y <- cbind(dd$y1, dd$y2, dd$y3)
M <- cbind(mu1, mu2, mu3)
S <- cbind(rep(s1, n), s2, rep(s3, n))
ll_dmvt <- vapply(seq_len(n), function(i) {
  D <- diag(S[i, ])
  mvtnorm::dmvt(Y[i, ], delta = M[i, ], sigma = D %*% C %*% D, df = nu,
                log = TRUE)
}, 0)
cat(sprintf("objective at shared point  %.12f\n", -nll_frm))
cat(sprintf("sum mvtnorm::dmvt          %.12f\n", sum(ll_dmvt)))
cat(sprintf("relative residual          %.3e\n",
            abs(-nll_frm - sum(ll_dmvt)) / abs(sum(ll_dmvt))))

# rescor_row_loglik() at the optimum reproduces logLik() row by row
dpv <- frmtmb:::eval_dpars(fit)
rl <- rescor_row_loglik(fit, dpv)
cat(sprintf("logLik(fit)                %.12f\n", as.numeric(logLik(fit))))
cat(sprintf("sum rescor_row_loglik      %.12f\n", sum(rl)))

## 2. independent RTMB objective, optimized separately
nll_ref <- function(q) {
  "c" <- RTMB::ADoverload("c")
  m1 <- q$b1[1] + q$b1[2] * x
  m2 <- q$b2[1] + q$b2[2] * x
  m3 <- q$b3[1] + q$b3[2] * x
  sd1 <- exp(q$ls1)
  sd2 <- exp(q$ls2[1] + q$ls2[2] * z)
  sd3 <- exp(q$ls3)
  Zs <- RTMB::matrix(c((dd$y1 - m1) / sd1, (dd$y2 - m2) / sd2,
                       (dd$y3 - m3) / sd3), n, K)
  Cq <- frmtmb::us_chol_cor(q$tr, K)
  nuq <- 1 + exp(q$lnu)
  Ci <- RTMB::solve(Cq)
  qv <- as.vector(((Zs %*% Ci) * Zs) %*% rep(1, K))
  ld <- -2 * (RTMB::dmvnorm(rep(0, K), 0, Cq, log = TRUE) +
                0.5 * K * log(2 * pi))
  ll <- n * (lgamma((nuq + K) / 2) - lgamma(nuq / 2) -
               K / 2 * log(nuq * pi) - ld / 2) -
    (nuq + K) / 2 * sum(log(1 + qv / nuq)) -
    n * log(sd1) - sum(log(sd2)) - n * log(sd3)
  -ll
}
ob2 <- RTMB::MakeADFun(nll_ref, list(b1 = c(0, 0), b2 = c(0, 0),
                                     b3 = c(0, 0), ls1 = 0, ls2 = c(0, 0),
                                     ls3 = 0, tr = numeric(3), lnu = 1),
                       silent = TRUE)
op2 <- nlminb(ob2$par, ob2$fn, ob2$gr,
              control = list(iter.max = 2000, eval.max = 2000))
cat(sprintf("ML logLik frmtmb           %.10f\n", as.numeric(logLik(fit))))
cat(sprintf("ML logLik reference        %.10f\n", -op2$objective))
nu_ref <- 1 + exp(op2$par[["lnu"]])
nu_frm <- 1 + exp(fit$estimates$betad[match("nu_(Intercept)", dn)])
cat(sprintf("nu frmtmb %.6f  reference %.6f\n", nu_frm, nu_ref))
cat(sprintf("rescor frmtmb: %s\n",
            paste(sprintf("%.6f", rescor_matrix(fit)[upper.tri(C)]),
                  collapse = " ")))
Cr <- frmtmb::us_chol_cor(op2$par[names(op2$par) == "tr"], K)
cat(sprintf("rescor reference: %s\n",
            paste(sprintf("%.6f", Cr[upper.tri(Cr)]), collapse = " ")))

## 3. large nu: the student objective approaches the gaussian one
fg <- frm(bf(y1 ~ x) + bf(y2 ~ x, sigma ~ z) + bf(y3 ~ x) +
            set_rescor(TRUE) + gaussian(), data = dd)
pg <- fg$opt$par
ps <- c(pg[names(pg) == "beta"], NA, pg[names(pg) == "thetar"])
# the student template has nu_(Intercept) inside betad; build it by name
bd_s <- numeric(length(dn))
bd_g <- pg[names(pg) == "betad"]
dng <- names(fg$frame$par_template$betad)
bd_s[match(dng, dn)] <- bd_g
for (lnu in c(5, 10, 20)) {
  bd_s[match("nu_(Intercept)", dn)] <- lnu
  pp <- c(pg[names(pg) == "beta"], bd_s, pg[names(pg) == "thetar"])
  cat(sprintf("log(nu - 1) = %2d: student %.10f gaussian %.10f\n", lnu,
              -obj$fn(pp), -fg$obj$fn(pg)))
}

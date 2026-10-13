# Reviewer of lane seflags, re-check: where does degree 7's SE/lmer of
# 1.09 come from? The exact beta Hessian of the ML LMM objective at
# lmer()'s estimates, X' V^-1 X with V = sigma^2 I + tau^2 Z Z', its
# unit-diagonal spectrum against the lane's re-measured one, and its
# inverse against lmer()'s vcov(), and the lane's reported SEs against
# both.
#   Rscript dev/seflags-rev2-deg7exact.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev2-common.R")
set.seed(7)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
  stats::rnorm(200, 0, 0.2)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
rhs <- paste0("x", 1:7, collapse = " + ")
fo <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
m <- suppressWarnings(lme4::lmer(fo, data = d, REML = FALSE,
  control = lme4::lmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 1e5))))
X <- lme4::getME(m, "X")
Z <- as.matrix(lme4::getME(m, "Z"))
sig <- stats::sigma(m)
tau <- attr(lme4::VarCorr(m)$g, "stddev")[[1]]
V <- sig^2 * diag(nrow(X)) + tau^2 * Z %*% t(Z)
Hb <- crossprod(X, solve(V, X))
D <- sqrt(diag(Hb))
es <- eigen(Hb / outer(D, D), symmetric = TRUE)
ev <- es$values
invH <- function(H) { dd <- sqrt(diag(H)); e2 <- eigen(H / outer(dd, dd), symmetric = TRUE); (e2$vectors %*% (t(e2$vectors) / e2$values)) / outer(dd, dd) }
cat("exact unit-diagonal spectrum (absolute), smallest 3:", signif(rev(ev)[1:3], 4), "| largest", signif(max(ev), 4), "| smallest 3 / largest:",
    signif(rev(ev)[1:3] / max(ev), 4), "\n")
se_exact <- sqrt(diag(invH(Hb)))
ref <- sqrt(diag(as.matrix(stats::vcov(m))))
cat("exact (X'V^-1X)^-1 SE / lmer vcov SE:", format(se_exact / ref, digits = 6),
    "\n")
f <- suppressWarnings(frm(fo, data = d))
cn <- c("(Intercept)", paste0("x", 1:7))
se <- se_outer(f)[cn]
cat("lane SE / exact at lmer's point:", format(se / se_exact, digits = 4), "\n")
cat("lane sigma", exp(f$opt$par[grepl("^sigma", ns$outer_par_names(f))]),
    "lmer sigma", sig, "| lane sd_g", exp(f$opt$par[ns$outer_par_names(f) ==
                                                      "theta_1"]),
    "lmer", tau, "\n")
# the exact Hessian at frmtmb's own estimates
on <- ns$outer_par_names(f)
s2 <- exp(f$opt$par[grepl("^sigma", on)])
t2 <- exp(f$opt$par[on == "theta_1"])
V2 <- s2^2 * diag(nrow(X)) + t2^2 * Z %*% t(Z)
Hb2 <- crossprod(X, solve(V2, X))
cat("lane SE / exact at frmtmb's point:",
    format(se / sqrt(diag(invH(Hb2))), digits = 4), "\n")
# accurate to the square root of the condition number: QR of the
# whitened design, V = L'L, Xw = L'^-1 X, cov = (R'R)^-1
qr_se <- function(s, t) {
  Vv <- s^2 * diag(nrow(X)) + t^2 * Z %*% t(Z)
  L <- chol(Vv)
  Xw <- backsolve(L, X, transpose = TRUE)
  q <- qr(Xw, tol = 1e-20)
  R <- qr.R(q)
  Ri <- backsolve(R, diag(ncol(X)))
  v <- numeric(ncol(X))
  v[q$pivot] <- sqrt(rowSums(Ri^2))
  v
}
se_qr <- qr_se(sig, tau)
cat("QR (accurate) SE / lmer vcov SE:", format(se_qr / ref, digits = 6), "\n")
cat("lane SE / QR at lmer's point:", format(se / se_qr, digits = 5), "\n")
cat("lane SE / QR at frmtmb's point:", format(se / qr_se(s2, t2), digits = 5),
    "\n")
cat("lm SE / QR(lm) check:", format(sqrt(diag(vcov(stats::lm(
  stats::as.formula(paste("y ~", rhs)), data = d)))) /
  qr_se(stats::sigma(stats::lm(stats::as.formula(paste("y ~", rhs)),
                               data = d)), 0), digits = 6), "\n")

# Smoke test of sar() in frmtmb.spatial: lag and error, gaussian and
# student, logLik against a dense multivariate normal at the estimates.
# Seed 2, 7 x 7 rook lattice, row-standardized W.
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.spatial"),
    "\n")
set.seed(2)
grid <- expand.grid(r = 1:7, c = 1:7)
n <- nrow(grid)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(grid$r[i] - grid$r[j]) + abs(grid$c[i] - grid$c[j]) == 1)
})
W <- A / rowSums(A)
d <- data.frame(x = rnorm(n))
IA <- diag(n) - 0.5 * W
d$ylag <- as.vector(solve(IA, 1 + 0.5 * d$x + rnorm(n, 0, 0.7)))
d$yerr <- 1 + 0.5 * d$x + as.vector(solve(IA, rnorm(n, 0, 0.7)))
dense <- function(y, m, S) {
  R <- chol(S)
  z <- backsolve(R, y - m, transpose = TRUE)
  -0.5 * length(y) * log(2 * pi) - sum(log(diag(R))) - 0.5 * sum(z^2)
}
fl <- frm(ylag ~ x + sar(W), data = d, data2 = list(W = W))
print(summary(fl))
rho <- confint_varcorr(fl)
print(rho)
eta <- frm_linpred(fl)
s <- sigma(fl)
rl <- variables(fl)
cat("vars:", rl, "\n")
r_hat <- hypothesis(fl, "lagsar = 0", class = NULL)
print(r_hat)
lagsar <- unname(frmtmb:::autocor_natural(fl$estimates$thetaac,
                                          fl$frame$autocor[[1]]))
Ar <- diag(n) - lagsar * W
S <- s^2 * solve(crossprod(Ar))
cat("lag logLik", format(as.numeric(logLik(fl)), digits = 12), "dense",
    format(dense(d$ylag, solve(Ar, eta), S), digits = 12), "\n")
cat("fitted = A^-1 eta:", isTRUE(all.equal(unname(fitted(fl)[, 1]),
                                           as.vector(solve(Ar, eta)))),
    "\n")
fe <- frm(yerr ~ x + sar(W, type = "error"), data = d, data2 = list(W = W))
errorsar <- unname(frmtmb:::autocor_natural(fe$estimates$thetaac,
                                            fe$frame$autocor[[1]]))
Ae <- diag(n) - errorsar * W
cat("error logLik", format(as.numeric(logLik(fe)), digits = 12), "dense",
    format(dense(d$yerr, frm_linpred(fe), sigma(fe)^2 *
                   solve(crossprod(Ae))), digits = 12), "\n")
fs <- frm(ylag ~ x + sar(W), family = student(), data = d,
          data2 = list(W = W))
print(summary(fs))
print(head(residuals(fl, type = "pearson")))
print(head(simulate(fl, nsim = 2, seed = 1)))
print(dim(autocor_matrix(fl)))
r <- try(fitted(fl, newdata = d[1:3, ]))
r <- try(frm(bf(ylag ~ x + sar(W), sigma ~ x), data = d,
             data2 = list(W = W)))
r <- try(frm(ylag ~ x + sar(W, type = "mixed"), data = d,
             data2 = list(W = W)))
r <- try(frm(ylag ~ x + sar(W), family = poisson(), data = d,
             data2 = list(W = W)))
print(frm_compat("sar()", "gaussian"))
print(frm_compat("sar()", "poisson"))

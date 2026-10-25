# Reviewer: independent dense check of the sar() and fcor() log densities
# at points AWAY from the estimates, against mvtnorm::dmvnorm()/dmvt()
# with the implied mean and covariance. Gaussian and student at nu = 5.
# Rscript dev/review-sarfcor-dense.R
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb.spatial); library(mvtnorm)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(31)
k <- 6
g <- expand.grid(r = 1:k, c = 1:k)
n <- nrow(g)
Adj <- outer(1:n, 1:n, function(i, j)
  as.numeric(pmax(abs(g$r[i] - g$r[j]), abs(g$c[i] - g$c[j])) == 1))
W <- Adj / rowSums(Adj)  # row-standardized queen: complex rounding?
ev <- eigen(W, only.values = TRUE)$values
cat("W eigen complex:", is.complex(ev), " range:",
    format(range(Re(ev)), digits = 5), "\n")
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * d$x + rnorm(n)))
X <- cbind(1, d$x)
V <- 0.5^abs(outer(1:n, 1:n, "-")) + diag(0.2, n)

report <- function(lab, ours, ref) {
  cat(sprintf("%-26s frmtmb %.10f dense %.10f diff %.3g\n", lab, ours, ref,
              ours - ref))
}
setp <- function(fit, beta, betad, rho = NULL) {
  p <- fit$obj$env$last.par.best
  p[names(p) == "beta"] <- beta
  p[names(p) == "betad"] <- betad
  if (!is.null(rho)) {
    aux <- fit$frame[["autocor"]][[1L]][["aux"]]
    p[names(p) == "thetaac"] <- qlogis((rho - aux$lo) / (aux$hi - aux$lo))
  }
  p
}
b <- c(0.7, -0.3)
mu <- as.vector(X %*% b)
for (type in c("lag", "error")) {
  for (rho in c(-0.9, 0.3, 0.85)) {
    A <- diag(n) - rho * W
    m <- if (type == "lag") solve(A, mu) else mu
    S0 <- solve(crossprod(A))
    fit <- frm(bf(y ~ x + sar(W, type = type)), data = d,
               data2 = list(W = W))
    nm <- names(fit$obj$env$last.par.best)
    sg <- 1.3
    p <- setp(fit, b, log(sg), rho)
    report(paste("gauss", type, rho), -fit$obj$fn(p),
           dmvnorm(d$y, m, sg^2 * S0, log = TRUE))
    fs <- frm(bf(y ~ x + sar(W, type = type)), family = student(), data = d,
              data2 = list(W = W))
    for (nu in c(3, 5, 30)) {
      p <- setp(fs, b, c(log(sg), log(nu - 1)), rho)
      report(paste("student", type, rho, "nu", nu), -fs$obj$fn(p),
             dmvt(d$y, m, sg^2 * S0, df = nu, log = TRUE))
    }
  }
}
cat("sar lag fit names:", unique(nm), "\n")

# fcor: constant sigma, het sigma ~ z, student hom and het
fit <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
report("fcor gauss", -fit$obj$fn(setp(fit, b, log(0.8))),
       dmvnorm(d$y, mu, 0.8^2 * V, log = TRUE))
fh <- frm(bf(y ~ x + fcor(V), sigma ~ z), data = d, data2 = list(V = V))
bd <- c(-0.2, 0.3)
D <- diag(exp(bd[1] + bd[2] * d$z))
report("fcor het gauss", -fh$obj$fn(setp(fh, b, bd)),
       dmvnorm(d$y, mu, D %*% V %*% D, log = TRUE))
fs <- frm(y ~ x + fcor(V), family = student(), data = d, data2 = list(V = V))
for (nu in c(3, 5, 30)) {
  report(paste("fcor student nu", nu),
         -fs$obj$fn(setp(fs, b, c(log(0.8), log(nu - 1)))),
         dmvt(d$y, mu, 0.8^2 * V, df = nu, log = TRUE))
}
fsh <- frm(bf(y ~ x + fcor(V), sigma ~ z), family = student(), data = d,
           data2 = list(V = V))
report("fcor het student nu 5",
       -fsh$obj$fn(setp(fsh, b, c(bd, log(5 - 1)))),
       dmvt(d$y, mu, D %*% V %*% D, df = 5, log = TRUE))
# a predicted nu: brms writes code that does not type-check; frmtmb?
r <- tryCatch({
  f <- frm(bf(y ~ x + fcor(V), nu ~ z), family = student(), data = d,
           data2 = list(V = V))
  paste("FIT", format(logLik(f)))
}, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("fcor + nu ~ z:", r, "\n")
r <- tryCatch({
  f <- frm(bf(y ~ x + sar(W), nu ~ z), family = student(), data = d,
           data2 = list(W = W))
  paste("FIT", format(logLik(f)))
}, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat("sar + nu ~ z:", r, "\n")

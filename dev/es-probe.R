# esicar lane probe: is the new esicar the exact hard-constrained model?
lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
library(frmtmb)

car_data <- function() {
  w <- matrix(0, 16, 16)
  g <- expand.grid(r = 1:4, c = 1:4)
  for (i in 1:16) for (j in 1:16) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) w[i, j] <- 1
  }
  dimnames(w) <- list(paste0("L", 1:16), paste0("L", 1:16))
  set.seed(42)
  kmat <- diag(rowSums(w)) - w + matrix(1 / (1e-3 * 16)^2, 16, 16)
  phi <- 1.2 * drop(crossprod(chol(solve(kmat)), rnorm(16)))
  loc <- factor(rep(rownames(w), each = 6), levels = rownames(w))
  d <- data.frame(loc = loc, x = rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + rnorm(nrow(d), 0, 0.5)
  list(d = d, W = w)
}
s <- car_data(); d <- s$d; W <- s$W

fit_of <- function(ty, con_sd = 1e-3) {
  frm(bf(y ~ x + car(W, gr = loc, type = ty, con_sd = con_sd)) + gaussian(),
      data = d, data2 = list(W = W))
}
f_es <- fit_of("esicar"); f_ic <- fit_of("icar")

cat("== log-likelihoods ==\n")
cat(sprintf("esicar %.15f  df %d\n", as.numeric(logLik(f_es)),
            attr(logLik(f_es), "df")))
cat(sprintf("icar   %.15f  df %d\n", as.numeric(logLik(f_ic)),
            attr(logLik(f_ic), "df")))
cat(sprintf("gap (esicar - icar) %.6e\n",
            as.numeric(logLik(f_es)) - as.numeric(logLik(f_ic))))
cat(sprintf("theta esicar %.15f  icar %.15f\n",
            f_es$estimates$theta[[1]], f_ic$estimates$theta[[1]]))
cat(sprintf("sdcar esicar %.10f  icar %.10f  rel %.3e\n",
            exp(f_es$estimates$theta[[1]]), exp(f_ic$estimates$theta[[1]]),
            exp(f_es$estimates$theta[[1]]) / exp(f_ic$estimates$theta[[1]]) - 1))

cat("\n== sum-to-zero ==\n")
fld <- function(f) {
  r <- ranef(f)
  r <- if (is.list(r) && !is.data.frame(r[[1]]) && !is.matrix(r[[1]])) r[[1]] else r
  as.numeric(r[[1]][, 1])
}
fld_es <- fld(f_es); fld_ic <- fld(f_ic)
cat(sprintf("esicar sum(field) %.6e   max|field| %.6f\n",
            sum(fld_es), max(abs(fld_es))))
cat(sprintf("icar   sum(field) %.6e   max|field| %.6f\n",
            sum(fld_ic), max(abs(fld_ic))))
cat(sprintf("esicar raw b sum  %.6e\n", sum(f_es$estimates$b)))

cat("\n== con_sd invariance ==\n")
for (cs in c(1e-2, 1e-3, 1e-4, 1e-5)) {
  a <- as.numeric(logLik(fit_of("esicar", cs)))
  b <- as.numeric(logLik(fit_of("icar", cs)))
  cat(sprintf("con_sd %.0e  esicar %.12f  icar %.12f\n", cs, a, b))
}

cat("\n== independent hard-constrained reference ==\n")
# y = X beta + Z f + eps, f the ICAR field under an EXACT sum-to-zero
# constraint: cov(f) = sdcar^2 * pinv(L). Marginal ML, no frmtmb.
L <- diag(rowSums(W)) - W
ev <- eigen(L, symmetric = TRUE)
pos <- ev$values > 1e-8 * max(ev$values)
Lp <- ev$vectors[, pos] %*% diag(1 / ev$values[pos]) %*% t(ev$vectors[, pos])
X <- model.matrix(~ x, d)
Z <- model.matrix(~ loc - 1, d)
ZLZ <- Z %*% Lp %*% t(Z)
y <- d$y; n <- length(y)
nll <- function(p) {
  V <- exp(2 * p[1]) * ZLZ + exp(2 * p[2]) * diag(n)
  ch <- chol(V)
  Vi <- chol2inv(ch)
  bh <- solve(crossprod(X, Vi %*% X), crossprod(X, Vi %*% y))
  r <- y - X %*% bh
  0.5 * (n * log(2 * pi) + 2 * sum(log(diag(ch))) +
           sum(r * (Vi %*% r)))
}
o <- optim(c(log(1.4), log(0.5)), nll, method = "BFGS",
           control = list(reltol = 1e-14))
cat(sprintf("reference logLik %.15f  sdcar %.10f  sigma %.10f\n",
            -o$value, exp(o$par[1]), exp(o$par[2])))
cat(sprintf("esicar - reference   %.6e\n",
            as.numeric(logLik(f_es)) + o$value))
cat(sprintf("icar   - reference   %.6e\n",
            as.numeric(logLik(f_ic)) + o$value))
cat(sprintf("sdcar esicar/reference - 1  %.3e\n",
            exp(f_es$estimates$theta[[1]]) / exp(o$par[1]) - 1))

cat("\n== predict / simulate / VarCorr ==\n")
print(VarCorr(f_es))
set.seed(7)
sm <- simulate(f_es, nsim = 3)
cat("simulate dim:", paste(dim(sm), collapse = " x "), "\n")
p1 <- predict(f_es)
p2 <- predict(f_es, newdata = d[1:5, ])
cat(sprintf("predict head %s\n", paste(round(head(p1, 3), 6), collapse = " ")))
cat(sprintf("predict newdata %s\n", paste(round(p2, 6), collapse = " ")))
cat(sprintf("logLik(escar) %.10f  logLik(bym2) %.10f\n",
            as.numeric(logLik(fit_of("escar"))),
            as.numeric(logLik(fit_of("bym2")))))

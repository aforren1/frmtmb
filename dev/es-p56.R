# P5 (singleton component) and P6 (non-gaussian con_sd invariance).
lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
library(frmtmb)
marginal_ml <- function(y, X, Z, Sigma_fn, start) {
  nll <- function(p) {
    if (any(!is.finite(p)) || any(abs(p) > 25)) return(1e10)
    V <- Z %*% Sigma_fn(p) %*% t(Z) + exp(2 * p[length(p)]) * diag(length(y))
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    Xs <- backsolve(R, X, transpose = TRUE)
    bh <- solve(crossprod(Xs), crossprod(Xs, backsolve(R, y, transpose = TRUE)))
    r <- y - X %*% bh
    0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) +
             sum(backsolve(R, r, transpose = TRUE)^2))
  }
  o <- stats::optim(start, nll, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 800))
  list(logLik = -o$value, par = o$par)
}

# --- P5: two 4-node paths plus one isolated node -----------------------
n <- 9
W <- matrix(0, n, n)
for (i in c(1, 2, 3, 5, 6, 7)) { W[i, i + 1] <- 1; W[i + 1, i] <- 1 }
lv <- paste0("L", seq_len(n)); dimnames(W) <- list(lv, lv)
set.seed(5)
L <- diag(rowSums(W)) - W
ev <- eigen(L, symmetric = TRUE)
pos <- ev$values > 1e-8 * max(ev$values)
cat("rank(L):", sum(pos), " expected n - c =", n - 3, "\n")
Lp <- ev$vectors[, pos] %*% diag(1 / ev$values[pos]) %*% t(ev$vectors[, pos])
phi <- drop(crossprod(chol(Lp + diag(1e-10, n)), rnorm(n)))
loc <- factor(rep(lv, each = 8), levels = lv)
d <- data.frame(loc = loc, x = rnorm(length(loc)))
d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + rnorm(nrow(d), 0, 0.4)
fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar")) + gaussian(),
           data = d, data2 = list(W = W))
bk <- fit$frame$re_blocks[[1]]
cat("n_comp:", bk$aux_car$n_comp, " nj:", paste(bk$aux_car$nj), "\n")
re <- ranef(fit)[[1]][, 1]
cat(sprintf("comp sums: %.3e %.3e  singleton value: %s (identical to 0: %s)\n",
            sum(re[1:4]), sum(re[5:8]), format(re[[9]]),
            identical(re[[9]], 0)))
ref <- marginal_ml(d$y, model.matrix(~x, d), model.matrix(~ loc - 1, d),
                   function(p) exp(2 * p[1]) * Lp, c(0, log(0.4)))
cat(sprintf("esicar %.12f  reference %.12f  gap %.3e\n",
            as.numeric(logLik(fit)), ref$logLik,
            as.numeric(logLik(fit)) - ref$logLik))
fi <- frm(bf(y ~ x + car(W, gr = loc, type = "icar")) + gaussian(),
          data = d, data2 = list(W = W))
cat(sprintf("icar %.12f  gap to ref %.3e  its singleton %.3e\n",
            as.numeric(logLik(fi)), as.numeric(logLik(fi)) - ref$logLik,
            ranef(fi)[[1]][, 1][[9]]))
cat("escar on the same W: ")
cat(tryCatch({frm(bf(y ~ x + car(W, gr = loc, type = "escar")) + gaussian(),
                  data = d, data2 = list(W = W)); "ACCEPTED"},
             error = function(e) paste("refused:", substr(conditionMessage(e), 1, 60))),
    "\n")

# --- P6: poisson -------------------------------------------------------
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c)); m <- nrow(g)
  M <- matrix(0, m, m)
  for (i in seq_len(m)) for (j in seq_len(m)) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) M[i, j] <- 1
  }
  dimnames(M) <- list(paste0("L", seq_len(m)), paste0("L", seq_len(m)))
  M
}
set.seed(21)
W2 <- lattice_W(4, 4); n2 <- nrow(W2)
K2 <- diag(rowSums(W2)) - W2 + matrix(1 / (1e-3 * n2)^2, n2, n2)
ph <- 0.6 * drop(crossprod(chol(solve(K2)), rnorm(n2)))
lc <- factor(rep(rownames(W2), each = 6), levels = rownames(W2))
dp <- data.frame(loc = lc, x = rnorm(length(lc)))
dp$y <- rpois(nrow(dp), exp(1 + 0.3 * dp$x + ph[as.integer(dp$loc)]))
llp <- function(ty, cs) {
  as.numeric(logLik(frm(bf(y ~ x + car(W2, gr = loc, type = ty,
                                       con_sd = cs)) + poisson(),
                        data = dp, data2 = list(W2 = W2))))
}
cs <- c(1e-2, 1e-3, 1e-4)
es <- vapply(cs, function(z) llp("esicar", z), 0)
ic <- vapply(cs, function(z) llp("icar", z), 0)
for (i in seq_along(cs)) {
  cat(sprintf("poisson %.0e  esicar %.12f  icar %.12f\n", cs[i], es[i], ic[i]))
}
cat(sprintf("esicar spread %.3e   icar spread %.3e\n",
            diff(range(es)), diff(range(ic))))
fp <- frm(bf(y ~ x + car(W2, gr = loc, type = "esicar")) + poisson(),
          data = dp, data2 = list(W2 = W2))
cat(sprintf("poisson esicar field sum %.3e\n", sum(ranef(fp)[[1]])))

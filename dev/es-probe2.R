lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
library(frmtmb)
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c)); n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
  }
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}
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
                    control = list(reltol = 1e-14, maxit = 500))
  list(logLik = -o$value, par = o$par)
}
car_lattice_data <- function(seed, r = 4, c = 4, per = 6, sd_car = 1.2,
                             sigma = 0.5, con_sd = 1e-3) {
  set.seed(seed); W <- lattice_W(r, c); n <- nrow(W)
  L <- diag(rowSums(W)) - W
  K <- L + matrix(1 / (con_sd * n)^2, n, n)
  phi <- sd_car * drop(crossprod(chol(solve(K)), stats::rnorm(n)))
  loc <- factor(rep(rownames(W), each = per), levels = rownames(W))
  d <- data.frame(loc = loc, x = stats::rnorm(length(loc)))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + stats::rnorm(nrow(d), 0, sigma)
  list(d = d, W = W, L = L, K = K, n = n,
       X = stats::model.matrix(~x, d), Z = stats::model.matrix(~ loc - 1, d))
}
s <- car_lattice_data(42); W <- s$W
A <- rbind(diag(s$n - 1), -1)
Vz <- solve(t(A) %*% s$L %*% A)
hard <- marginal_ml(s$d$y, s$X, s$Z,
                    function(p) exp(2 * p[1]) * (A %*% Vz %*% t(A)), c(0, 0))
f <- function(ty, cs = 1e-3) frm(bf(y ~ x + car(W, gr = loc, type = ty,
                                                con_sd = cs)) + gaussian(),
                                 data = s$d)
fe <- f("esicar"); fi <- f("icar")
cat(sprintf("hard   %.15f  par %s\n", hard$logLik,
            paste(round(hard$par, 10), collapse = " ")))
cat(sprintf("esicar %.15f  theta %.12f\n", as.numeric(logLik(fe)),
            fe$estimates$theta[[1]]))
cat(sprintf("icar   %.15f  theta %.12f\n", as.numeric(logLik(fi)),
            fi$estimates$theta[[1]]))
cat(sprintf("esicar - hard %.4e   icar - hard %.4e   esicar - icar %.4e\n",
            as.numeric(logLik(fe)) - hard$logLik,
            as.numeric(logLik(fi)) - hard$logLik,
            as.numeric(logLik(fe)) - as.numeric(logLik(fi))))
cat(sprintf("esicar sum(ranef) %.3e   icar sum(ranef) %.3e\n",
            sum(ranef(fe)[[1]]), sum(ranef(fi)[[1]])))
cat("betad/theta esicar:", paste(round(c(fe$estimates$theta,
                                         fe$estimates$betad), 8),
                                 collapse = " "), "\n")
cat("con_sd sweep (esicar, icar):\n")
for (cs in c(1e-2, 1e-3, 1e-4, 1e-5)) {
  cat(sprintf("  %.0e  %.12f  %.12f\n", cs,
              as.numeric(logLik(f("esicar", cs))),
              as.numeric(logLik(f("icar", cs)))))
}
cat("\nSE of sdcar (esicar):\n")
print(head(confint(fe), 8))

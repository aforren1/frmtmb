# esicar alongside an ordinary (1 | g): expand_b() must center one block
# and copy the other, so the marginal likelihood is the sum of a
# hard-constrained ICAR and an iid intercept.
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
set.seed(11)
W <- lattice_W(4, 4); n <- nrow(W)
L <- diag(rowSums(W)) - W
ev <- eigen(L, symmetric = TRUE); pos <- ev$values > 1e-8 * max(ev$values)
Lp <- ev$vectors[, pos] %*% diag(1 / ev$values[pos]) %*% t(ev$vectors[, pos])
phi <- 1.0 * drop(crossprod(chol(Lp + diag(1e-10, n)), rnorm(n)))
phi <- phi - mean(phi)
loc <- factor(rep(rownames(W), each = 8), levels = rownames(W))
g2 <- factor(rep(paste0("g", 1:8), length.out = length(loc)))
u <- rnorm(8, 0, 0.6)
d <- data.frame(loc = loc, g2 = g2, x = rnorm(length(loc)))
d$y <- 1 + 0.5 * d$x + phi[as.integer(d$loc)] + u[as.integer(d$g2)] +
  rnorm(nrow(d), 0, 0.5)

fit <- frm(bf(y ~ x + car(W, gr = loc, type = "esicar") + (1 | g2)) +
             gaussian(), data = d, data2 = list(W = W))

X <- model.matrix(~ x, d)
Zl <- model.matrix(~ loc - 1, d)
Zg <- model.matrix(~ g2 - 1, d)
ZLZ <- Zl %*% Lp %*% t(Zl)
ZGZ <- Zg %*% t(Zg)
y <- d$y; N <- length(y)
nll <- function(p) {
  V <- exp(2 * p[1]) * ZLZ + exp(2 * p[2]) * ZGZ + exp(2 * p[3]) * diag(N)
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(1e10)
  Xs <- backsolve(R, X, transpose = TRUE)
  bh <- solve(crossprod(Xs), crossprod(Xs, backsolve(R, y, transpose = TRUE)))
  r <- y - X %*% bh
  0.5 * (N * log(2 * pi) + 2 * sum(log(diag(R))) +
           sum(backsolve(R, r, transpose = TRUE)^2))
}
o <- optim(c(0, log(0.6), log(0.5)), nll, method = "BFGS",
           control = list(reltol = 1e-14, maxit = 800))
cat(sprintf("frmtmb  %.12f\n", as.numeric(logLik(fit))))
cat(sprintf("reference %.12f\n", -o$value))
cat(sprintf("gap %.3e\n", as.numeric(logLik(fit)) + o$value))
cat("theta frmtmb:", paste(round(fit$estimates$theta, 8), collapse = " "), "\n")
cat("reference    :", paste(round(o$par[1:2], 8), collapse = " "), "\n")
re <- ranef(fit)
cat("blocks:", paste(names(re), collapse = " | "), "\n")
for (nm in names(re)) {
  v <- re[[nm]]
  cat(sprintf("  %s: n=%d sum=%.3e\n", nm, nrow(v), sum(v)))
}

## Shared data + objective definitions for the RTMB vs RTMBp comparison.
## Sourced by every benchmark process. Deliberately free of any RTMB/RTMBp
## reference: the objective closes over sparse matrices built here, and the
## generics it calls (dnorm/dpois/%*%/sum) resolve against whichever of the
## two packages the caller attached. That is the "identical code modulo the
## namespace" requirement.

## ---- shape A: InstEval LMM, y ~ service + (1|s) + (1|d) --------------
build_A <- function() {
  d <- lme4::InstEval
  X  <- Matrix::sparse.model.matrix(~ service, d)
  Z1 <- Matrix::t(Matrix::fac2sparse(d$s))
  Z2 <- Matrix::t(Matrix::fac2sparse(d$d))
  list(
    name = "A_InstEval_LMM",
    y = as.numeric(d$y), X = X, Z1 = Z1, Z2 = Z2,
    n = nrow(X), q1 = ncol(Z1), q2 = ncol(Z2),
    parameters = list(beta = c(mean(as.numeric(d$y)), 0),
                      logsigma = 0, logs1 = -1, logs2 = -1,
                      u1 = numeric(ncol(Z1)), u2 = numeric(ncol(Z2))),
    random = c("u1", "u2"))
}

make_f_A <- function(dat) {
  y <- dat$y; X <- dat$X; Z1 <- dat$Z1; Z2 <- dat$Z2
  function(p) {
    beta <- p$beta; u1 <- p$u1; u2 <- p$u2
    sigma <- exp(p$logsigma); s1 <- exp(p$logs1); s2 <- exp(p$logs2)
    eta <- as.vector(X %*% beta) + as.vector(Z1 %*% u1) + as.vector(Z2 %*% u2)
    -sum(dnorm(y, eta, sigma, TRUE)) -
      sum(dnorm(u1, 0, s1, TRUE)) - sum(dnorm(u2, 0, s2, TRUE))
  }
}

## ---- shape B: Poisson GLMM, n = 3e5, 1000 groups ---------------------
## Accumulation-heavy: the inner Hessian for a single (1|g) is diagonal
## (1000 x 1000), so the sparse Cholesky is nearly free and the
## per-observation sum dominates. This is where parallel accumulation
## should show up if it shows up anywhere.
build_B <- function() {
  set.seed(20260901)
  n <- 300000L; ng <- 1000L
  g <- factor(sample.int(ng, n, replace = TRUE), levels = seq_len(ng))
  x <- rnorm(n)
  u <- rnorm(ng, 0, 0.5)
  eta <- 0.5 + 0.3 * x + u[as.integer(g)]
  y <- rpois(n, exp(eta))
  X <- Matrix::sparse.model.matrix(~ x, data.frame(x = x))
  Z <- Matrix::t(Matrix::fac2sparse(g))
  list(
    name = "B_Poisson_GLMM",
    y = y, X = X, Z = Z, n = n, ng = ng,
    parameters = list(beta = c(0, 0), logsu = -1, u = numeric(ng)),
    random = "u")
}

make_f_B <- function(dat) {
  y <- dat$y; X <- dat$X; Z <- dat$Z
  function(p) {
    beta <- p$beta; u <- p$u
    su <- exp(p$logsu)
    eta <- as.vector(X %*% beta) + as.vector(Z %*% u)
    -sum(dpois(y, exp(eta), TRUE)) - sum(dnorm(u, 0, su, TRUE))
  }
}

build_shape <- function(shape) if (shape == "A") build_A() else build_B()
make_f      <- function(shape, dat) if (shape == "A") make_f_A(dat) else make_f_B(dat)

tm <- function(expr) {
  t <- proc.time()[["elapsed"]]; force(expr)
  proc.time()[["elapsed"]] - t
}

## Elapsed plus process CPU time. On Windows proc.time() sums CPU over all
## threads of the process, so cpu/elapsed > 1 is direct evidence that
## OpenMP threads actually ran. This is what separates "threads
## configured" from "threads used".
tmc <- function(expr) {
  t <- proc.time(); force(expr); d <- proc.time() - t
  c(elapsed = d[["elapsed"]], cpu = d[["user.self"]] + d[["sys.self"]])
}

# Item 3: whose defect is it? The lane's hand tape, then: (1) no atomic product;
# (2) no atomic Cholesky either (a Cholesky by hand on the tape); (3) a plain TMB
# C++ template with no RTMB; (4) other inner orders. No frmtmb. Seed 20261102.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
library(RTMB)
cat("lib:", dirname(system.file(package = "RTMB")), "TMB", as.character(packageVersion("TMB")),
    "RTMB", as.character(packageVersion("RTMB")), "\n")
set.seed(20261102)
n0 <- 60L
x0 <- round(runif(n0, 0, 6), 3); e0 <- rnorm(n0)
y <- (sin(2 * x0) + 1e-2 * e0)[1:40]; x <- x0[1:40]; n <- 40L
u <- sort(unique(x)); m <- length(u); J <- match(x, u); D2 <- outer(u, u, "-")^2
handchol <- function(K) {
  # a right-looking Cholesky on advectors, one op at a time, no atomic
  "[<-" <- RTMB::ADoverload("[<-")
  mm <- nrow(K)
  L <- K * 0
  for (j in seq_len(mm)) {
    s <- K[j, j]
    if (j > 1L) for (k in seq_len(j - 1L)) s <- s - L[j, k] * L[j, k]
    L[j, j] <- sqrt(s)
    if (j < mm) for (i in (j + 1L):mm) {
      s <- K[i, j]
      if (j > 1L) for (k in seq_len(j - 1L)) s <- s - L[i, k] * L[j, k]
      L[i, j] <- s / L[j, j]
    }
  }
  L
}
mk <- function(atomic_chol, atomic_mm) function(p) {
  K <- exp(-D2 / (2 * exp(2 * p$lrho))) + diag(1e-12, m)
  L <- if (atomic_chol) t(chol(K)) else handchol(K)
  if (atomic_mm) { o <- TapeConfig(); TapeConfig(matmul = "atomic") }
  v <- as.vector(L %*% p$b)
  if (atomic_mm) TapeConfig(o)
  f <- exp(p$lsd) * v
  eta <- p$mu + f[J]
  -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) - sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
}
at <- function(o) {
  p <- o$par
  p[names(p) == "lsd"] <- 1.08396; p[names(p) == "lrho"] <- 0.36733; p[names(p) == "ls"] <- -4.57404
  p
}
ng <- function(o, p, h = 1e-4) vapply(seq_along(p), function(i) {
  e <- replace(numeric(length(p)), i, h)
  (o$fn(p + e) - o$fn(p - e)) / (2 * h)
}, 0)
for (ac in c(TRUE, FALSE)) for (am in c(TRUE, FALSE)) {
  fn <- mk(ac, am)
  for (ord in c("mu_first", "b_first")) for (comp in 0:1) {
    old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
    TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
    pars <- if (ord == "mu_first") list(mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3, b = numeric(m)) else
      list(b = numeric(m), mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3)
    o <- MakeADFun(fn, pars, random = c("mu", "b"), silent = TRUE)
    force(o$env$spHess)
    TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB")
    p <- at(o)
    g <- o$gr(p)
    cat(sprintf("RTD atomicchol %d atomicmm %d %-8s compress %d fn %.8f gr lrho %.6f\n", ac, am, ord, comp, o$fn(p),
                g[names(p) == "lrho" | seq_along(g) == 3L][1]))
    flush.console()
  }
}
cat("== which orders of the inner vector, atomic chol and product, compress 1\n")
fn <- mk(TRUE, TRUE)
perm <- function(ix) {
  # the field's position order as data: permute u
  u2 <- u[ix]
  D2p <- outer(u2, u2, "-")^2
  Jp <- match(x, u2)
  function(p) {
    K <- exp(-D2p / (2 * exp(2 * p$lrho))) + diag(1e-12, m)
    L <- t(chol(K))
    o <- TapeConfig(); TapeConfig(matmul = "atomic")
    v <- as.vector(L %*% p$b); TapeConfig(o)
    f <- exp(p$lsd) * v
    eta <- p$mu + f[Jp]
    -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) - sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
  }
}
set.seed(77)
farthest <- local({
  d2 <- D2; p <- integer(m); p[1L] <- which.min(rowSums(d2)); gap <- d2[, p[1L]]; gap[p[1L]] <- -1
  for (k in seq_len(m - 1L) + 1L) { p[k] <- which.max(gap); gap <- pmin(gap, d2[, p[k]]); gap[p[seq_len(k)]] <- -1 }
  p
})
for (lab in c("sorted", "farthest", "random1", "random2", "reverse")) {
  ix <- switch(lab, sorted = seq_len(m), farthest = farthest, random1 = sample(m), random2 = sample(m), reverse = rev(seq_len(m)))
  fp <- perm(ix)
  for (ord in c("mu_first", "b_first")) {
    res <- list()
    for (comp in 0:1) {
      old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
      TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
      pars <- if (ord == "mu_first") list(mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3, b = numeric(m)) else
        list(b = numeric(m), mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3)
      o <- MakeADFun(fp, pars, random = c("mu", "b"), silent = TRUE)
      force(o$env$spHess)
      TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB")
      res[[comp + 1L]] <- o$gr(at(o))
    }
    cat(sprintf("RTD order %-9s %-8s plain %s | compressed %s | maxdiff %.2e\n", lab, ord,
                paste(sprintf("%.6f", res[[1]]), collapse = " "), paste(sprintf("%.6f", res[[2]]), collapse = " "),
                max(abs(res[[1]] - res[[2]]))))
    flush.console()
  }
}
cat("== plain TMB template, no RTMB\n")
tdir <- file.path(tempdir(), "rev"); dir.create(tdir)
cpp <- file.path(tdir, "revdefect.cpp")
writeLines(c(
  "#include <TMB.hpp>",
  "template<class Type> Type objective_function<Type>::operator() () {",
  "  DATA_VECTOR(y); DATA_IVECTOR(J); DATA_MATRIX(D2);",
  "  PARAMETER(mu); PARAMETER(ls); PARAMETER(lsd); PARAMETER(lrho); PARAMETER_VECTOR(b);",
  "  int m = D2.rows(); matrix<Type> K(m, m);",
  "  for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) K(i, j) = exp(-D2(i, j) / (Type(2) * exp(Type(2) * lrho))) + (i == j ? Type(1e-12) : Type(0));",
  "  Eigen::LLT< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > llt(K);",
  "  matrix<Type> L = llt.matrixL();",
  "  vector<Type> v = L * b.matrix();",
  "  Type nll = -dnorm(b, Type(0), Type(1), true).sum();",
  "  for (int i = 0; i < y.size(); i++) nll -= dnorm(y(i), mu + exp(lsd) * v(J(i)), exp(ls), true);",
  "  return nll;",
  "}"), cpp)
r <- tryCatch({ TMB::compile(cpp, framework = "TMBad"); TRUE }, error = function(e) { cat("compile failed:", conditionMessage(e), "\n"); FALSE })
if (r) {
  dyn.load(TMB::dynlib(sub("[.]cpp$", "", cpp)))
  for (comp in 0:1) {
    TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "revdefect")
    o <- TMB::MakeADFun(list(y = y, J = J - 1L, D2 = D2), list(mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3, b = numeric(m)),
                        random = c("mu", "b"), DLL = "revdefect", silent = TRUE)
    force(o$env$spHess)
    p <- at(o)
    g <- o$gr(p)
    fd <- ng(o, p)
    cat(sprintf("RTD TMB-C++ compress %d fn %.8f gr %s | fd %s\n", comp, o$fn(p), paste(sprintf("%.6f", g), collapse = " "),
                paste(sprintf("%.6f", fd), collapse = " ")))
  }
  TMB::config(tmbad.sparse_hessian_compress = 0L, DLL = "revdefect")
}
cat("DONE\n")

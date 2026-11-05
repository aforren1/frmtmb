# Walks the hand tape toward frm()'s to find what the compressed Hessian's reverse
# sweep gets wrong with the intercept in the inner vector. No frmtmb. Seed 20261102.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
library(RTMB)
set.seed(20261102)
n0 <- 60L
x0 <- round(runif(n0, 0, 6), 3)
e0 <- rnorm(n0)
y <- (sin(2 * x0) + 1e-2 * e0)[1:40]
x <- x0[1:40]
n <- 40L
u <- sort(unique(x))
m <- length(u)
J <- match(x, u)
D2 <- outer(u, u, "-")^2
Z <- Matrix::sparseMatrix(i = seq_len(n), j = J, x = 1, dims = c(n, m))
Z <- methods::as(Z, "CsparseMatrix")
X <- matrix(1, n, 1)
fld <- function(p) {
  K <- exp(-D2 / (2 * exp(2 * p$lrho))) + diag(1e-12, m)
  L <- t(chol(K))
  o <- TapeConfig(); TapeConfig(matmul = "atomic")
  v <- as.vector(L %*% p$b)
  TapeConfig(o)
  (exp(p$lsd) + 1e-150) * v
}
variants <- list(
  vector = function(p) {
    f <- fld(p)
    eta <- p$mu + f[J]
    -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) - sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
  },
  Z = function(p) {
    f <- fld(p)
    eta <- p$mu + as.vector(Z %*% f)
    -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) - sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
  },
  Xbeta = function(p) {
    f <- fld(p)
    eta <- as.vector(X %*% p$mu) + as.vector(Z %*% f)
    -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) - sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
  },
  cvec = function(p) {
    "[<-" <- RTMB::ADoverload("[<-")
    cv <- rep(p$b[1] * 0, m)
    cv[seq_len(m)] <- fld(p)
    eta <- as.vector(X %*% p$mu) + as.vector(Z %*% cv)
    -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) - sum(RTMB::dnorm(y, eta, exp(p$ls), log = TRUE))
  })
ng <- function(o, p, h = 1e-4) vapply(seq_along(p), function(i) {
  e <- replace(numeric(length(p)), i, h)
  (o$fn(p + e) - o$fn(p - e)) / (2 * h)
}, 0)
for (vn in names(variants)) for (comp in 0:1) {
  old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
  o <- MakeADFun(variants[[vn]], list(mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3, b = numeric(m)),
                 random = c("mu", "b"), silent = TRUE)
  force(o$env$spHess)
  TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB")
  p <- o$par; p[names(p) == "lsd"] <- 1.08396; p[names(p) == "lrho"] <- 0.36733; p[names(p) == "ls"] <- -4.57404
  g <- o$gr(p); fd <- ng(o, p)
  cat(sprintf("NCGPRB %-8s compress %d fn %.8f gr %s | fd %s | Hnodes %d\n", vn, comp, o$fn(p),
              paste(format(g, digits = 5), collapse = " "), paste(format(fd, digits = 5), collapse = " "),
              nrow(GetTape(o, "ADHess")$data.frame())))
}
cat("--- order and size\n")
for (ord in c("mu_first", "b_first")) for (comp in 0:1) {
  old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
  pars <- if (ord == "mu_first") list(mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3, b = numeric(m)) else
    list(ls = log(0.01), lsd = 0.5, lrho = 0.3, b = numeric(m), mu = 0)
  o <- MakeADFun(variants$vector, pars, random = c("mu", "b"), silent = TRUE)
  force(o$env$spHess)
  TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB")
  p <- o$par; p[names(p) == "lsd"] <- 1.08396; p[names(p) == "lrho"] <- 0.36733; p[names(p) == "ls"] <- -4.57404
  g <- o$gr(p); fd <- ng(o, p)
  cat(sprintf("NCGPRB %-8s compress %d fn %.8f gr %s | fd %s\n", ord, comp, o$fn(p),
              paste(format(g, digits = 5), collapse = " "), paste(format(fd, digits = 5), collapse = " ")))
}

# REML by hand on the tape: the intercept in the inner vector beside z. Outer
# gradient against a central difference, by compress x matmul, at low noise, where
# frm()'s tape read d/d lscale -0.34 against 0.07. No frmtmb. Seed 20261102.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
library(RTMB)
set.seed(20261102)
n0 <- 60L
x <- round(runif(n0, 0, 6), 3)[1:40]
y <- (sin(2 * round(x, 3)))
set.seed(20261102)
x <- round(runif(n0, 0, 6), 3)
e0 <- rnorm(n0)
y <- (sin(2 * x) + 1e-2 * e0)[1:40]
x <- x[1:40]
n <- 40L
D2 <- outer(x, x, "-")^2
fun <- function(p) {
  K <- exp(-D2 / (2 * exp(2 * p$lrho))) + RTMB::diag(1e-12, n)
  L <- t(chol(K))
  f <- exp(p$lsd) * as.vector(L %*% p$b)
  -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) -
    sum(RTMB::dnorm(y, p$mu + f, exp(p$ls), log = TRUE))
}
ng <- function(o, p, h = 1e-4) vapply(seq_along(p), function(i) {
  e <- replace(numeric(length(p)), i, h)
  (o$fn(p + e) - o$fn(p - e)) / (2 * h)
}, 0)
for (rand in list(c("b", "mu"), "b")) for (comp in 0:1) for (mm in c("plain", "atomic")) {
  old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  tc <- TapeConfig()
  TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "RTMB")
  TapeConfig(matmul = mm)
  o <- MakeADFun(fun, list(mu = 0, lsd = 0.5, lrho = 0.3, ls = log(0.01), b = numeric(n)),
                 random = rand, silent = TRUE)
  force(o$env$spHess)
  TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB")
  TapeConfig(tc)
  p <- o$par; p[names(p) == "lsd"] <- 1.08396; p[names(p) == "lrho"] <- 0.36733; p[names(p) == "ls"] <- -4.57404
  g <- o$gr(p); fd <- ng(o, p)
  cat(sprintf("NCGPRA random %-7s compress %d matmul %-6s fn %.8f max |gr - fd| %.2e  gr %s fd %s\n",
              paste(rand, collapse = ","), comp, mm, o$fn(p), max(abs(g - fd)),
              paste(format(g, digits = 4), collapse = " "), paste(format(fd, digits = 4), collapse = " ")))
}

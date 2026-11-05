# A failed factor under the tape settings frm() uses: atomic product, compressed
# Hessian, and frm()'s own tape with its nugget replaced. Which of them returns a
# number at a length scale where chol() fails? Seed 20261103.
#   Rscript dev/ncgp-cholfail2.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
n <- 60L
set.seed(20261102)
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
y <- sin(2 * x) + 0.1 * e0
u <- sort(unique(x)); J <- match(x, u); n <- length(u)
D2 <- outer(u, u, "-")^2
nug <- 1e-17
fails <- function(r) {
  K <- exp(-D2 / (2 * r^2)); diag(K) <- diag(K) + nug
  inherits(try(chol(K), silent = TRUE), "try-error")
}
cat("NCGPC2 chol at rho 5 fails:", fails(5), " rho 0.01 fails:", fails(0.01), "\n")
mk <- function(atomic, compress) {
  old <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
  on.exit(TMB::config(tmbad.sparse_hessian_compress = old, DLL = "RTMB"))
  TMB::config(tmbad.sparse_hessian_compress = as.integer(compress), DLL = "RTMB")
  MakeADFun(function(p) {
    r <- exp(p$lrho)
    K <- exp(-D2 / (2 * r^2)) + RTMB::diag(nug, n)
    L <- t(chol(K))
    if (atomic) {
      o <- TapeConfig(); TapeConfig(matmul = "atomic")
      f <- as.vector(L %*% p$b)
      TapeConfig(o)
    } else f <- as.vector(L %*% p$b)
    f <- exp(p$lsd) * f
    -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) -
      sum(RTMB::dnorm(y, p$mu + f[J], exp(p$ls), log = TRUE))
  }, list(mu = 0, lsd = 0, lrho = log(0.01), ls = 0, b = numeric(n)), random = "b", silent = TRUE)
}
for (a in c(FALSE, TRUE)) for (cc in c(FALSE, TRUE)) {
  obj <- mk(a, cc)
  pb <- c(0, 0, log(5), 0)
  v <- tryCatch(obj$fn(pb), error = function(e) conditionMessage(e))
  g <- tryCatch(obj$gr(pb), error = function(e) conditionMessage(e))
  # at the same point, the closed form at nugget 1e-12 for scale
  cat(sprintf("NCGPC2 atomic %d compress %d: fn %s gr %s\n", a, cc, format(v, digits = 8),
              paste(format(g, digits = 5), collapse = " ")))
  obj$fn(c(0, 0, log(1), 0))
  # the inner solution there
  pb2 <- obj$env$last.par
  cat("   last.par b head", head(pb2[names(pb2) == "b"], 3), "\n")
}
library(frmtmb)
d <- data.frame(x = x, y = y)
f2 <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
ns <- asNamespace("frmtmb")
unlockBinding("gp_nugget", ns); assign("gp_nugget", nug, envir = ns)
fr <- f2$frame
bk <- fr$re_blocks[[1]]
K2 <- frmtmb:::gp_corr(c(0, log(5)), bk)
cat("NCGPC2 frm gp_corr rho 5 fails:", inherits(try(chol(K2), silent = TRUE), "try-error"), "\n")
tp <- MakeTape(function(th) { K <- frmtmb:::gp_corr(th, bk); t(chol(K))[nrow(K), nrow(K)] }, c(0, log(0.01)))
cat("NCGPC2 is.numeric on an advector:", is.numeric(RTMB::advector(1)), "\n")
cat("NCGPC2 raw tape at rho 5:", tryCatch(tp(c(0, log(5))), error = function(e) paste("ERROR", conditionMessage(e))), "\n")
obj <- frmtmb:::gp_tape_scope(fr, MakeADFun(frmtmb:::build_objective(fr), f2$estimates, random = "b", silent = TRUE))
bad <- f2$opt$par; bad[which(names(bad) == "theta")[2]] <- log(5)
cat("NCGPC2 frm tape:", obj$fn(bad), "gr", obj$gr(bad), "\n")
print(bad)
# frm() tape from a start where the factor already fails when the tape is built
est2 <- f2$estimates; est2$theta[2] <- log(5)
r <- tryCatch({
  o2 <- frmtmb:::gp_tape_scope(fr, MakeADFun(frmtmb:::build_objective(fr), est2, random = "b", silent = TRUE))
  c(o2$fn(o2$par), o2$fn(f2$opt$par))
}, error = function(e) paste("ERROR", conditionMessage(e)))
cat("NCGPC2 tape built at a failing theta:", r, "\n")

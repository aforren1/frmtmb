# What does a failed Cholesky factor do on the AD tape? Hand closure, no
# frmtmb, seed 20261103. A relative 1e-12 never failed on the scan
# (dev/ncgp-cholscan.R), so the failure is constructed: nugget 1e-16, where the
# factor fails at most length scales. Reads obj$fn and obj$gr before and after,
# with RTMB's chol(), to see whether an error, a NaN or a wrong number comes out.
#   Rscript dev/ncgp-cholfail.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
library(RTMB)
n <- 200L
set.seed(20261103)
x <- sort(runif(n, 0, 10))
y <- sin(x) + 0.3 * rnorm(n)
D2 <- outer(x, x, "-")^2
nug <- 1e-16
fails <- function(r) {
  K <- exp(-D2 / (2 * r^2))
  diag(K) <- diag(K) + nug
  inherits(try(chol(K), silent = TRUE), "try-error")
}
grid <- exp(seq(log(0.05), log(100), length.out = 200))
res <- vapply(grid, fails, NA)
cat("NCGPC fails at", sum(res), "of", length(res), "; first", grid[which(res)[1]],
    "last ok", grid[max(which(!res))], "\n")
good <- grid[which(!res)[1]]
bad <- grid[which(res)[length(which(res)) %/% 2L]]
cat("NCGPC good rho", good, "bad rho", bad, "\n")
obj <- MakeADFun(function(p) {
  r <- exp(p$lrho)
  K <- exp(-D2 / (2 * r^2)) + RTMB::diag(nug, n)
  L <- t(chol(K))
  f <- exp(p$lsd) * as.vector(L %*% p$b)
  -sum(RTMB::dnorm(p$b, 0, 1, log = TRUE)) -
    sum(RTMB::dnorm(y, p$mu + f, exp(p$ls), log = TRUE))
}, list(mu = 0, lsd = 0, lrho = log(good), ls = 0, b = numeric(n)), random = "b",
silent = TRUE)
show <- function(tag, p) {
  vf <- tryCatch(obj$fn(p), error = function(e) paste("ERROR", conditionMessage(e)))
  vg <- tryCatch(obj$gr(p), error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("NCGPC %s: fn %s | gr %s\n", tag, paste(format(vf, digits = 10), collapse = " "),
              paste(format(vg, digits = 5), collapse = " ")))
}
pg <- c(0, 0, log(good), 0)
pb <- c(0, 0, log(bad), 0)
show("good", pg)
show("bad", pb)
show("good again", pg)
o <- tryCatch(nlminb(c(0, 0, log(bad) - 0.05, 0), obj$fn, obj$gr, control = list(iter.max = 50)),
              error = function(e) paste("ERROR", conditionMessage(e)))
if (is.list(o)) cat("NCGPC nlminb from near bad:", o$message, o$objective, exp(o$par[3]), "\n") else cat("NCGPC nlminb:", o, "\n")

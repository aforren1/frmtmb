# How often are the third derivatives of RTMB::pbeta() non-finite over the
# (q, a, b) an xbeta fit visits, for the direct form I_q(a, b) and for
# the complement 1 - I_{1-q}(b, a)? Seed 20260929, 4000 points,
# q = kappa / (1 + 2 kappa) with log kappa ~ U(-6, 2), log a and log b
# ~ U(-3, 6).
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
set.seed(20260929)
N <- 4000
lk <- runif(N, -6, 2); la <- runif(N, -3, 6); lb <- runif(N, -3, 6)
fd <- function(p) {
  k <- exp(p[1]); q <- k / (1 + 2 * k)
  log(pbeta(q, exp(p[2]), exp(p[3])))
}
fc <- function(p) {
  k <- exp(p[1]); q <- k / (1 + 2 * k)
  log1p(-pbeta(1 - q, exp(p[3]), exp(p[2])))
}
mk <- function(f) {
  F <- MakeTape(f, c(-1, 0, 0))
  J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
  list(F = F, J = J, H = H, T3 = T3)
}
D <- mk(fd); C <- mk(fc)
res <- t(vapply(seq_len(N), function(i) {
  x <- c(lk[i], la[i], lb[i])
  ref <- pbeta(exp(x[1]) / (1 + 2 * exp(x[1])), exp(x[2]), exp(x[3]),
               log.p = TRUE)
  one <- function(M) {
    c(val_ok = is.finite(M$F(x)) &&
        abs(M$F(x) - ref) <= 1e-10 * max(1, abs(ref)),
      g_ok = all(is.finite(M$J(x))), h_ok = all(is.finite(M$H(x))),
      t_ok = all(is.finite(M$T3(x))))
  }
  c(ref = ref, one(D), one(C))
}, numeric(9)))
colnames(res) <- c("ref", paste0("d_", c("val", "g", "h", "t")),
                   paste0("c_", c("val", "g", "h", "t")))
res <- as.data.frame(res)
cat("points", N, "\n")
cat("reference log P finite:", sum(is.finite(res$ref)), "\n")
fin <- is.finite(res$ref) & res$ref > -700
cat("of which log P > -700:", sum(fin), "\n")
for (nm in names(res)[-1]) {
  cat(sprintf("%-6s ok %4d of %4d\n", nm, sum(res[[nm]][fin] == 1), sum(fin)))
}
cat("third derivative finite in at least one form:",
    sum((res$d_t == 1 | res$c_t == 1)[fin]), "\n")
bad <- fin & res$d_t == 0
cat("direct-form NaN third derivative, by region:\n")
print(summary(data.frame(q = (exp(lk) / (1 + 2 * exp(lk)))[bad],
                         a = exp(la)[bad], b = exp(lb)[bad],
                         logP = res$ref[bad])))

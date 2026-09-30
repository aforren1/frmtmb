# The switched continued fraction for log I_x(a, b) on the RTMB tape:
# value against stats::pbeta(log.p = TRUE), gradient against a central
# difference of that reference, and whether second and third
# derivatives are finite, over x < 1/2 and shapes from e^-4 to e^16.
# Seed 20260929. Compared with RTMB::pbeta() on the same points.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(RTMB)
CFN <- as.integer(Sys.getenv("FAMS2_CFN", "50"))
log_ibeta <- function(x, a, b) frmtmb:::log_ibeta_half(x, a, b, N = CFN)
cat("continued fraction steps N =", CFN, "\n")
N_SWEEP <- as.integer(Sys.getenv("FAMS2_N", "3000"))
set.seed(20260929)
lk <- runif(N_SWEEP, -12, 3)
la <- runif(N_SWEEP, -4, 16)
lb <- runif(N_SWEEP, -4, 16)
# parameters: log kappa, log a, log b
ours <- function(p) {
  k <- exp(p[1]); x <- k / (1 + 2 * k)
  log_ibeta(x, exp(p[2]), exp(p[3]))
}
theirs <- function(p) {
  k <- exp(p[1]); x <- k / (1 + 2 * k)
  log(pbeta(x, exp(p[2]), exp(p[3])))
}
mk <- function(f) {
  F <- MakeTape(f, c(-1, 0, 0))
  J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
  list(F = F, J = J, H = H, T3 = T3)
}
O <- mk(ours); R <- mk(theirs)
refv <- function(p) {
  k <- exp(p[1])
  stats::pbeta(k / (1 + 2 * k), exp(p[2]), exp(p[3]), log.p = TRUE)
}
out <- t(vapply(seq_len(N_SWEEP), function(i) {
  p <- c(lk[i], la[i], lb[i])
  ref <- suppressWarnings(refv(p))
  h <- 1e-6
  gref <- vapply(1:3, function(j) {
    e <- replace(numeric(3), j, h)
    suppressWarnings((refv(p + e) - refv(p - e)) / (2 * h))
  }, 0)
  one <- function(M) {
    v <- M$F(p); g <- M$J(p)
    c(verr = abs(v - ref) / max(1, abs(ref)),
      gerr = max(abs(g - gref)) / max(1, abs(gref)),
      hfin = all(is.finite(M$H(p))), tfin = all(is.finite(M$T3(p))))
  }
  c(ref = ref, one(O), one(R))
}, numeric(9)))
colnames(out) <- c("ref", "o_verr", "o_gerr", "o_hfin", "o_tfin",
                   "r_verr", "r_gerr", "r_hfin", "r_tfin")
out <- as.data.frame(out)
keep <- is.finite(out$ref) & out$ref > -600
cat("points", N_SWEEP, "; with -600 < log P:", sum(keep), "\n")
s <- function(v) sprintf("max %.2e q99 %.2e q999 %.2e", max(v),
                         quantile(v, 0.99), quantile(v, 0.999))
cat("value rel err   ours  ", s(out$o_verr[keep]), "\n")
cat("value rel err   RTMB  ", s(out$r_verr[keep]), "\n")
cat("grad  rel err   ours  ", s(out$o_gerr[keep]), "(against a central difference)\n")
cat("grad  rel err   RTMB  ", s(out$r_gerr[keep]), "\n")
cat("hessian finite  ours", sum(out$o_hfin[keep]), " RTMB", sum(out$r_hfin[keep]), "\n")
cat("third   finite  ours", sum(out$o_tfin[keep]), " RTMB", sum(out$r_tfin[keep]), "\n")
bad <- keep & (out$o_verr > 1e-10 | out$o_tfin == 0)
if (any(bad)) {
  print(head(data.frame(x = exp(lk) / (1 + 2 * exp(lk)), a = exp(la),
                        b = exp(lb), out)[bad, ], 20))
}
# the second derivative against a central difference of our gradient
set.seed(7)
idx <- sample(which(keep), 200)
herr <- vapply(idx, function(i) {
  p <- c(lk[i], la[i], lb[i]); h <- 1e-5
  fd <- sapply(1:3, function(j) {
    e <- replace(numeric(3), j, h)
    (O$J(p + e) - O$J(p - e)) / (2 * h)
  })
  max(abs(O$H(p) - fd)) / max(1, abs(fd))
}, 0)
cat("hessian vs central difference of the gradient, 200 points:", s(herr), "\n")
# value error by the size of the larger shape
mx <- pmax(la, lb)[keep]
bin <- cut(exp(mx), c(0, 1e2, 1e4, 1e6, Inf))
cat("value rel err by max(a, b):\n")
print(tapply(out$o_verr[keep], bin, function(v) sprintf("max %.1e (%d)", max(v), length(v))))
print(tapply(out$r_verr[keep], bin, function(v) sprintf("RTMB max %.1e", max(v))))
xx <- exp(lk) / (1 + 2 * exp(lk)); aa <- exp(la); bb <- exp(lb)
side <- ifelse(xx < (aa + 1) / (aa + bb + 2), "direct", "complement")
o <- order(-out$o_verr * keep)[1:12]
print(data.frame(x = xx, a = aa, b = bb, side = side, m = (aa + 1) / (aa + bb + 2),
                 ref = out$ref, err = out$o_verr)[o, ])
cat("worst by side:\n")
print(tapply(out$o_verr[keep], side[keep], max))
# ours against RTMB's own AD derivatives, which are exact where finite
dd <- t(vapply(which(keep), function(i) {
  p <- c(lk[i], la[i], lb[i])
  gR <- R$J(p); hR <- R$H(p)
  c(g = max(abs(O$J(p) - gR)) / max(1, abs(gR)),
    h = if (all(is.finite(hR))) max(abs(O$H(p) - hR)) / max(1, abs(hR)) else NA)
}, numeric(2)))
cat("gradient vs RTMB::pbeta AD:", s(dd[, "g"]), "\n")
cat("hessian  vs RTMB::pbeta AD:", s(na.omit(dd[, "h"])), "\n")
cat("gradient vs RTMB AD by max(a, b):\n")
print(tapply(dd[, "g"], bin, function(v) sprintf("max %.1e", max(v))))

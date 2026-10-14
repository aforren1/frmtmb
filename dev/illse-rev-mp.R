# Reviewer of lane illse: shared setup and 256-bit helpers, written
`%||%` <- function(a, b) if (is.null(a)) b else a
# independently of the lane's dev/illse-mp.R.
#   REV_ARM = lane (default) | base
rev_arm <- Sys.getenv("REV_ARM", "lane")
.libPaths(c(if (rev_arm == "lane") "C:/Users/adf44/source/r/wt-illse-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(Rmpfr)
})
ns <- asNamespace("frmtmb")
cat(sprintf("arm %s | frmtmb from %s | R %s\n", rev_arm,
            find.package("frmtmb"), R.home()))
PREC <- 256L
m <- function(x) mpfr(x, PREC)
mM <- function(A) mpfrArray(as.numeric(A), PREC, dim = dim(as.matrix(A)))
num <- function(x) asNumeric(x)
# Gauss-Jordan inverse with partial pivoting in 256 bits
mp_inv <- function(A) {
  p <- nrow(A)
  B <- mM(diag(p))
  for (k in seq_len(p)) {
    piv <- k - 1L + which.max(num(abs(A[k:p, k])))
    if (piv != k) {
      t1 <- A[k, ]; A[k, ] <- A[piv, ]; A[piv, ] <- t1
      t2 <- B[k, ]; B[k, ] <- B[piv, ]; B[piv, ] <- t2
    }
    s <- A[k, k]
    A[k, ] <- A[k, ] / s
    B[k, ] <- B[k, ] / s
    for (i in seq_len(p)[-k]) {
      f <- A[i, k]
      A[i, ] <- A[i, ] - f * A[k, ]
      B[i, ] <- B[i, ] - f * B[k, ]
    }
  }
  B
}
# central second differences in 256 bits; step 1e-20 leaves truncation
# near 1e-40 and rounding near 1e-35 relative
mp_hessian <- function(f, x, h = 1e-20) {
  p <- length(x)
  x <- m(x)
  h <- m(h)
  H <- mM(matrix(0, p, p))
  ei <- function(i) { v <- m(rep(0, p)); v[i] <- h; v }
  f0 <- f(x)
  for (i in seq_len(p)) {
    H[i, i] <- (f(x + ei(i)) - 2 * f0 + f(x - ei(i))) / h^2
    if (i < p) for (j in (i + 1L):p) {
      v <- (f(x + ei(i) + ei(j)) - f(x + ei(i) - ei(j)) -
              f(x - ei(i) + ei(j)) + f(x - ei(i) - ei(j))) / (4 * h^2)
      H[i, j] <- v
      H[j, i] <- v
    }
  }
  H
}
# gradient of a scalar f by central differences in 256 bits
mp_grad <- function(f, x, h = 1e-30) {
  p <- length(x)
  x <- m(x)
  h <- m(h)
  out <- m(rep(0, p))
  for (i in seq_len(p)) {
    e <- m(rep(0, p)); e[i] <- h
    out[i] <- (f(x + e) - f(x - e)) / (2 * h)
  }
  out
}
# a' V a for each row of A (double), V mpfr
mp_rowvar <- function(A, V) {
  A <- as.matrix(A)
  Am <- mM(A)
  num(rowSums((Am %*% V) * Am))
}
# every entry of V (double) against Ve (mpfr): largest error over
# sqrt(Ve_ii Ve_jj), largest relative error over all entries, and over
# the diagonal
cov_err <- function(V, Ve) {
  V <- as.matrix(V)
  if (!all(is.finite(V))) return(c(cor = NaN, rel = NaN, diag = NaN))
  D <- mM(V) - Ve
  dg <- sqrt(num(diag(Ve)))
  c(cor = max(abs(num(D)) / outer(dg, dg)),
    rel = max(abs(num(D / Ve))),
    diag = max(abs(num(diag(D) / diag(Ve)))))
}
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "rev_err")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
wtag <- function(w) {
  if (!length(w)) return("-")
  k <- c(prod = "product such as", getv = "marginaleffects forms",
         lost = "not available", ill = "ill-conditioned",
         conv = "did not report", grad = "gradient", flat = "flat")
  hit <- names(k)[vapply(k, function(s) any(grepl(s, w, fixed = TRUE)), NA)]
  if (!length(hit)) "other" else paste(hit, collapse = ",")
}
ratio_str <- function(r) {
  if (!length(r)) return("NA")
  if (all(is.na(r))) return(if (all(is.nan(r))) "NaN" else "NA")
  sprintf("%.6f-%.6f (max dev %.1e)", min(r, na.rm = TRUE), max(r, na.rm = TRUE),
          max(abs(r - 1), na.rm = TRUE))
}

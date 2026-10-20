# Review of lane adhess: copy of illse2-common.R (lane illse2) with REV_ARM = lane|base|merged|illse2.
#   ILLSE2_ARM = lane (default) puts the lane library first;
#   base measures rellib-r9 (frmtmb 0.71.0) alone.
`%||%` <- function(a, b) if (is.null(a)) b else a
arm <- Sys.getenv("REV_ARM", "lane")
lib_lane <- c(lane = "C:/Users/adf44/source/r/wt-adhess-lib", merged = "C:/Users/adf44/source/r/adhess-rev-lib", illse2 = "C:/Users/adf44/source/r/wt-illse2-lib")[arm]
lib_base <- "C:/Users/adf44/source/r/rellib-r9"
lib_user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (arm == "r8") lib_base <- "C:/Users/adf44/source/r/rellib-r8"
.libPaths(c(if (arm != "base") lib_lane, lib_base, lib_user))
suppressMessages({
  library(frmtmb)
  library(Rmpfr)
})
ns <- asNamespace("frmtmb")
cat(sprintf("arm %s | frmtmb %s from %s | BLAS %s\n", arm,
            utils::packageVersion("frmtmb"),
            dirname(find.package("frmtmb")),
            basename(extSoftVersion()[["BLAS"]] %||% "")))
PREC <- 256L
m <- function(x) mpfr(x, PREC)
mM <- function(A) mpfrArray(as.numeric(A), PREC, dim = dim(as.matrix(A)))
num <- function(x) asNumeric(x)
col1 <- function(v) mpfr2array(v, dim = c(length(v), 1L))
# Gauss-Jordan inverse with partial pivoting in 256 bits
mp_inv <- function(A) {
  if (!inherits(A, "mpfrMatrix")) A <- mM(A)
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
# q' A q for an mpfr matrix A and a double vector q
mp_quad <- function(A, q) {
  qq <- col1(m(q))
  num(sum(qq * (A %*% qq)))
}
# negative log-likelihood of a gaussian fit in (beta, log sigma), exact
# in 256 bits through sufficient statistics
gauss_nll <- function(X, y) {
  Xm <- mM(X); ym <- m(y)
  XX <- crossprod(Xm); Xy <- crossprod(Xm, ym); yy <- sum(ym^2)
  p <- ncol(X); n <- length(y)
  function(th) {
    b <- th[seq_len(p)]
    s2 <- exp(2 * th[p + 1L])
    bm <- col1(b)
    rss <- yy - 2 * sum(Xy * b) + sum(bm * (XX %*% bm))
    (n * log(s2) + rss / s2) / 2
  }
}
# the same with (1 | g): beta, log sigma, log sd in that order
lmm_nll <- function(X, y, g) {
  Xm <- mM(X); ym <- m(y)
  XX <- crossprod(Xm); Xy <- crossprod(Xm, ym); yy <- sum(ym^2)
  p <- ncol(X)
  gi <- split(seq_along(y), g)
  S <- lapply(gi, function(ix) list(s = colSums(Xm[ix, , drop = FALSE]),
                                     sy = sum(ym[ix]), n = length(ix)))
  function(th) {
    b <- th[seq_len(p)]; a2 <- exp(2 * th[p + 1L]); c2 <- exp(2 * th[p + 2L])
    bm <- col1(b)
    rss <- yy - 2 * sum(Xy * b) + sum(bm * (XX %*% bm))
    out <- rss / a2
    for (st in S) {
      sr <- st$sy - sum(st$s * b)
      out <- out + (st$n - 1) * log(a2) + log(a2 + st$n * c2) -
        c2 * sr^2 / (a2 * (a2 + st$n * c2))
    }
    out / 2
  }
}
# exact Hessian of a gaussian nll in closed form (no second differences)
gauss_hess <- function(X, y, th) {
  p <- ncol(X)
  Xm <- mM(X)
  r <- m(y) - Xm %*% col1(m(th[seq_len(p)]))
  s2 <- exp(2 * m(th[p + 1L]))
  H <- mM(matrix(0, p + 1L, p + 1L))
  H[seq_len(p), seq_len(p)] <- crossprod(Xm) / s2
  hb <- 2 * crossprod(Xm, r) / s2
  H[seq_len(p), p + 1L] <- hb
  H[p + 1L, seq_len(p)] <- hb
  H[p + 1L, p + 1L] <- 2 * sum(r^2) / s2
  H
}
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "i2_err")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
is_err <- function(x) inherits(x, "i2_err")
wtag <- function(w) {
  if (!length(w)) return("-")
  k <- c(prod = "product such as", getv = "marginaleffects forms",
         lost = "not available", ill = "ill-conditioned",
         conv = "did not report", grad = "gradient", flat = "flat",
         short = "short", bound = "boundary")
  hit <- names(k)[vapply(k, function(s) any(grepl(s, w, fixed = TRUE)), NA)]
  if (!length(hit)) "other" else paste(hit, collapse = ",")
}
ratio_str <- function(r) {
  if (!length(r)) return("NA")
  if (all(is.na(r))) return(if (all(is.nan(r))) "NaN" else "NA")
  sprintf("%.9f-%.9f (max dev %.1e)", min(r, na.rm = TRUE),
          max(r, na.rm = TRUE), max(abs(r - 1), na.rm = TRUE))
}
lost_of <- function(f) {
  l <- ns$sdr_of(f)$se_lost
  if (!length(l)) return("none")
  tb <- table(l)
  paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = ", ")
}

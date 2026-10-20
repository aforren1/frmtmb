# Lane adhess: shared setup and exact references.
#   ADHESS_ARM = lane (default) puts the lane library first;
#   base measures rellib-r9 (frmtmb 0.71.0) alone.
`%||%` <- function(a, b) if (is.null(a)) b else a
arm <- Sys.getenv("ADHESS_ARM", "lane")
lib_lane <- "C:/Users/adf44/source/r/wt-adhess-lib"
lib_base <- "C:/Users/adf44/source/r/rellib-r9"
lib_user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (arm == "lane") lib_lane, lib_base, lib_user))
suppressMessages({
  library(frmtmb)
  library(Rmpfr)
})
ns <- asNamespace("frmtmb")
cat(sprintf("arm %s | frmtmb %s from %s | BLAS %s | LAPACK %s\n", arm,
            utils::packageVersion("frmtmb"),
            dirname(find.package("frmtmb")),
            basename(extSoftVersion()[["BLAS"]] %||% ""),
            basename(La_library() %||% "")))
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
# Exact information X' W X of a canonical-link GLM at the coefficients
# b, in 256 bits; the covariance is its inverse, rounded once to double
glm_info_exact <- function(X, b, family = c("poisson", "bernoulli")) {
  family <- match.arg(family)
  Xm <- mM(X)
  eta <- Xm %*% col1(m(b))
  w <- if (family == "poisson") exp(eta) else {
    mu <- 1 / (1 + exp(-eta))
    mu * (1 - mu)
  }
  Xw <- Xm * as.vector(w)
  crossprod(Xm, Xw)
}
glm_cov_exact <- function(X, b, family) {
  num(mp_inv(glm_info_exact(X, b, family)))
}
# the error measures of the round: SE ratio and the largest covariance
# entry error on the correlation scale, both against Vr
cov_err <- function(V, Vr) {
  dg <- sqrt(diag(Vr))
  c(se = max(abs(sqrt(diag(V)) / dg - 1)),
    cor = max(abs(V - Vr) / outer(dg, dg)))
}
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "ah_err")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(x) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
is_err <- function(x) inherits(x, "ah_err")
wtag <- function(w) {
  if (!length(w)) return("-")
  k <- c(prod = "product such as", getv = "marginaleffects forms",
         lost = "not available", ill = "ill-conditioned",
         conv = "did not report", grad = "gradient", flat = "flat",
         short = "short", bound = "boundary")
  hit <- names(k)[vapply(k, function(s) any(grepl(s, w, fixed = TRUE)), NA)]
  if (!length(hit)) "other" else paste(hit, collapse = ",")
}
lost_of <- function(f) {
  l <- ns$sdr_of(f)$se_lost
  if (!length(l)) return("none")
  tb <- table(l)
  paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = ", ")
}

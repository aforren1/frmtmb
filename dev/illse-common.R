# Lane illse: shared setup. ILLSE_ARM = "lane" (default) puts the lane
# library first; "base" measures rellib-r8 (frmtmb 0.70.0) alone; "dev"
# loads the worktree itself over rellib-r8 (pkgload::load_all()).
arm <- Sys.getenv("ILLSE_ARM", "lane")
lib_lane <- "C:/Users/adf44/source/r/wt-illse-lib"
lib_base <- "C:/Users/adf44/source/r/rellib-r8"
lib_user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (arm == "lane") lib_lane, lib_base, lib_user))
if (arm == "dev") {
  # the worktree as it is, without an install (pure R)
  suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-illse",
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
cat(sprintf("arm %s | frmtmb %s from %s | BLAS %s\n", arm,
            utils::packageVersion("frmtmb"),
            dirname(find.package("frmtmb")),
            basename(extSoftVersion()[["BLAS"]] %||% "")))
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-mp.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
# value and every warning of an expression
cap <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "illse_err")),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  list(v = v, w = w)
}
lost_of <- function(f) {
  l <- f$sdr$se_lost %||% frmtmb:::sdr_of(f)$se_lost
  if (!length(l)) return("none")
  tb <- table(l)
  paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = ", ")
}
# exact covariance of a gaussian fit without random effects at the
# fit's own point: the Hessian of n log s + RSS / (2 s^2) in
# (beta, log s), in 256-bit arithmetic, inverted there
gauss_exact_cov <- function(X, y, beta, logs) {
  p <- ncol(X)
  Xm <- mpfrArray(X, mp_prec, dim = dim(X))
  r <- mp(y) - Xm %*% mp(beta)
  s2 <- exp(2 * mp(logs))
  H <- mpfrArray(0, mp_prec, dim = c(p + 1, p + 1))
  H[1:p, 1:p] <- crossprod(Xm) / s2
  hb <- 2 * crossprod(Xm, r) / s2
  H[1:p, p + 1] <- hb
  H[p + 1, 1:p] <- hb
  H[p + 1, p + 1] <- 2 * sum(r^2) / s2
  mp_solve(H)
}
# exact covariance of a gaussian (1 | g) fit at the fit's own point: the
# marginal negative log-likelihood (which the Laplace objective is,
# exactly, for a gaussian model) in (beta, log sigma, log sd), with its
# Hessian by second differences in 256-bit arithmetic (mp_hess())
lmm_exact_cov <- function(X, y, g, par) {
  p <- ncol(X)
  Xm <- mpfrArray(X, mp_prec, dim = dim(X))
  ym <- mp(y)
  gi <- split(seq_along(y), g)
  # sufficient statistics, exact in 256 bits: the residual sum of
  # squares and each group's residual sum are polynomials in beta
  XX <- crossprod(Xm)
  Xy <- crossprod(Xm, ym)
  yy <- sum(ym^2)
  S <- lapply(gi, function(ix) {
    list(s = colSums(Xm[ix, , drop = FALSE]), sy = sum(ym[ix]),
         n = length(ix))
  })
  nll <- function(th) {
    b <- th[1:p]
    s2 <- exp(2 * th[p + 1])
    t2 <- exp(2 * th[p + 2])
    bm <- mpfr2array(b, dim = c(p, 1L))
    rss <- yy - 2 * sum(Xy * b) + sum(bm * (XX %*% bm))
    out <- (length(y) - length(gi)) * log(s2) + rss / s2
    for (st in S) {
      den <- s2 + st$n * t2
      rs <- st$sy - sum(st$s * b)
      out <- out + log(den) - t2 / (den * s2) * rs^2
    }
    out / 2
  }
  mp_solve(mp_hess(nll, par))
}

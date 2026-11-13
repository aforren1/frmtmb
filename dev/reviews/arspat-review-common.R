# Review of lane arspat: shared setup. Source from the worktree root.
lib <- Sys.getenv("ARSPAT_LIB", "C:/Users/adf44/source/r/wt-arspat-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
ns <- asNamespace("frmtmb")

chk <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", conditionMessage(w), "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) {
           cat("MSG[", paste(class(m)[1], collapse = ","), "]: ",
               conditionMessage(m), sep = "")
           invokeRestart("muffleMessage")
         }),
       error = function(e) { cat("ERROR: ", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) print(r)
  invisible(r)
}

# natural parameters of the ar() family term of a fit
ac_nat <- function(fit) {
  ac <- fit$frame$autocor[[1L]]
  th <- fit$estimates$thetaac[ac$theta_idx]
  ns$autocor_natural(th, ac)
}

# autocovariance of an ARMA with innovation sd sigma at lags 0..L
arma_gamma <- function(ar, ma, sigma, L) {
  ar <- if (length(ar)) ar else numeric(0)
  ma <- if (length(ma)) ma else numeric(0)
  psi <- c(1, stats::ARMAtoMA(ar, ma, 5000))
  g0 <- sigma^2 * sum(psi^2)
  rho <- stats::ARMAacf(ar, ma, lag.max = L)
  g0 * as.numeric(rho)
}

# the dense marginal log-density: every group's covariance over its time
# indices (all rows, missing included), observed rows kept
dense_margin <- function(y, mu, sig, tix, grp, ar, ma, se = NULL,
                         nu = NULL, cosy = NULL) {
  tot <- 0
  for (i in split(seq_along(y), grp)) {
    tt <- tix[i]
    L <- max(tt) - min(tt)
    if (!is.null(cosy)) {
      R <- matrix(cosy, length(i), length(i)); diag(R) <- 1
      S <- R * outer(sig[i], sig[i])
    } else {
      gm <- arma_gamma(ar, ma, 1, max(L, 1))
      S <- matrix(gm[abs(outer(tt, tt, "-")) + 1L], length(i)) *
        outer(sig[i], sig[i])
    }
    if (!is.null(se)) S <- S + diag(se[i]^2, length(i))
    ok <- !is.na(y[i])
    if (!any(ok)) next
    S <- S[ok, ok, drop = FALSE]
    tot <- tot + if (is.null(nu)) {
      mvtnorm::dmvnorm(y[i][ok], mu[i][ok], S, log = TRUE)
    } else {
      mvtnorm::dmvt(y[i][ok], mu[i][ok], S, df = nu, log = TRUE)
    }
  }
  tot
}

# Review of lane adhess: shared setup.
#   REV_ARM = lane | base | merged | illse2
`%||%` <- function(a, b) if (is.null(a)) b else a
arm <- Sys.getenv("REV_ARM", "lane")
libs <- list(
  lane = "C:/Users/adf44/source/r/wt-adhess-lib",
  merged = "C:/Users/adf44/source/r/adhess-rev-lib",
  illse2 = "C:/Users/adf44/source/r/wt-illse2-lib",
  base = NULL)
.libPaths(c(libs[[arm]], "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat(sprintf("arm %s | frmtmb %s from %s\n", arm,
            utils::packageVersion("frmtmb"), dirname(find.package("frmtmb"))))
cap <- function(expr) {
  w <- character()
  m <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "rev_err")),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(x) {
      m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")
    })
  list(v = v, w = w, m = m)
}
# Richardson Jacobian of the AD gradient at `at`, objective state kept
rich_hess <- function(obj, at, d = 1e-3, r = 6) {
  st <- ns$obj_state_save(obj)
  on.exit(ns$obj_state_restore(obj, st))
  J <- numDeriv::jacobian(function(p) as.numeric(obj$gr(p)), at,
                          method = "Richardson",
                          method.args = list(eps = 1e-4, d = d,
                                             zero.tol = 1e-6, r = r, v = 2))
  (J + t(J)) / 2
}
cov_err <- function(V, Vr) {
  dg <- sqrt(diag(Vr))
  ok <- is.finite(diag(V)) & is.finite(dg) & dg > 0
  if (!any(ok)) return(c(se = NA, cor = NA))
  V <- V[ok, ok, drop = FALSE]
  Vr <- Vr[ok, ok, drop = FALSE]
  dg <- dg[ok]
  c(se = max(abs(sqrt(diag(V)) / dg - 1)),
    cor = max(abs(V - Vr) / outer(dg, dg)))
}
the_point <- function(fit) {
  u <- fit$par_units
  scaled <- !(is.null(u) || all(u == 1))
  if (scaled) fit$opt$par else ns$sdr_outer_point(fit$obj)
}

# Reviewer of lane setier, re-check: dev/setier-rev2-window.R with the
# probe replaced (this process only) by symmetric second differences:
# keep a direction when c(full) = f(+d) + f(-d) - 2 f0 >= 2 grad_tol and
# c(full) / c(half) is in (3, 5.5). The linear term a fit short of exact
# stationarity leaves cancels in c().
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
sym <- function(fit, p, free, dir_free, lambda) {
  tol <- fit$control$grad_tol %||% 1e-3
  obj <- fit$obj
  saved <- obj_state_save(obj)
  on.exit(obj_state_restore(obj, saved), add = TRUE)
  d <- numeric(length(p))
  d[free] <- dir_free * sqrt(4 * tol / lambda)
  ev <- function(x) tryCatch(obj$fn(x), error = function(e) NA_real_)
  f0 <- ev(p)
  cf <- ev(p + d) + ev(p - d) - 2 * f0
  ch <- ev(p + d / 2) + ev(p - d / 2) - 2 * f0
  is.finite(cf) && is.finite(ch) && cf >= 2 * tol && ch > 0 &&
    cf / ch > se_quad_lo && cf / ch < se_quad_hi
}
environment(sym) <- ns
unlockBinding("se_curvature_real", ns)
assign("se_curvature_real", sym, envir = ns)
commandArgs <- function(...) c("sym", "20")
source("dev/setier-rev2-window.R")

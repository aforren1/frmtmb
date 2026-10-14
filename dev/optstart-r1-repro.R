# Lane optstart, item R1: a copy of lane remeasure's
# dev/remeasure-repro-missings.R (that lane's findings, R1), unchanged
# below this header.
# Two outcomes of the brms_missings vignette that changed at 0.70.0 with
# the optimum unchanged: fit_imp2 now reports false convergence, and
# fit_imp3 now withholds the standard errors of both residual SDs.
#
#   Rscript dev/remeasure-repro-missings.R <lib>
#
# <lib> is the library to measure (the lane's 0.67.0 build, or
# rellib-r8). The data are the vignette's: mice::nhanes, and for
# fit_imp3 a measurement-error column drawn under the port's own seed for
# that expression (port_seed("brms_missings.10.1"), dev/brms-port/shim.R).
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
source("dev/brms-port/shim.R")
data("nhanes", package = "mice")
w <- character()
keep <- function(expr) {
  withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
}
show <- function(fit, label) {
  cat("\n##", label, "\n")
  cat(sprintf("logLik %.9f  convergence %d  message %s\n",
              as.numeric(logLik(fit)), fit$opt$convergence,
              fit$opt$message %||% ""))
  cat("warnings:", if (length(w)) paste(substr(w, 1, 160), collapse = " | ")
      else "none", "\n")
  se <- suppressWarnings(sqrt(diag(vcov(fit))))
  cat("vcov() standard errors:\n")
  print(signif(se, 5))
  p <- fit$opt$par
  g <- fit$obj$gr(p)
  cat(sprintf("max |gradient| at the optimum %.3g\n", max(abs(g))))
  # An independent reference: the observed information of the objective
  # the optimizer minimized, by finite differences of its gradient
  # (stats::optimHess), inverted. frmtmb's SE machinery is not involved.
  H <- stats::optimHess(p, fit$obj$fn, fit$obj$gr)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  cat(sprintf("finite-difference Hessian: %d parameters, eigenvalues %s\n",
              length(p), paste(signif(range(ev), 4), collapse = " to ")))
  ref <- tryCatch(sqrt(diag(solve(H))), error = function(e) NA)
  print(signif(stats::setNames(ref, names(p)), 5))
  cat("diagnose():\n")
  print(tryCatch(diagnose(fit), error = function(e) conditionMessage(e)))
  if (!is.null(fit$cache$se_explained)) {
    cat("fit$cache$se_explained:\n"); str(fit$cache$se_explained)
  }
  w <<- character()
}
`%||%` <- function(a, b) if (is.null(a)) b else a
bform2 <- bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
  set_rescor(FALSE)
f2 <- keep(frm(bform2, data = nhanes))
show(f2, "fit_imp2: bmi | mi() ~ age * mi(chl), chl | mi() ~ age")
set.seed(port_seed("brms_missings.10.1"))
nhanes$se <- rexp(nrow(nhanes), 2)
bform3 <- bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi(se) ~ age) +
  set_rescor(FALSE)
f3 <- keep(frm(bform3, data = nhanes))
show(f3, "fit_imp3: the same with chl | mi(se)")
# Which 0.70.0 change: the same two fits with the per-parameter unit
# scaling turned off, where the build has the switch.
if ("autoscale" %in% names(formals(frmtmb_control))) {
  ctl <- frmtmb_control(autoscale = FALSE)
  show(keep(frm(bform2, data = nhanes, control = ctl)),
       "fit_imp2, autoscale = FALSE")
  show(keep(frm(bform3, data = nhanes, control = ctl)),
       "fit_imp3, autoscale = FALSE")
}

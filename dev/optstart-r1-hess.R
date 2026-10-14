# Lane optstart, item R1: which standard errors of brms_missings'
# fit_imp2 and fit_imp3 are right. At the fit's optimum: vcov() (the
# fit's own route, the outer Hessian in natural units), the inverse of
# stats::optimHess() in raw units (0.67.0's route, and lane remeasure's
# reference), and the inverse of a Richardson-extrapolated Jacobian of
# the AD gradient (numDeriv), the most accurate of the three.
#   Rscript dev/optstart-r1-hess.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
source("dev/brms-port/shim.R")
data("nhanes", package = "mice")
set.seed(port_seed("brms_missings.10.1"))
se_col <- rexp(nrow(nhanes), 2)
forms <- list(
  fit_imp2 = bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
    set_rescor(FALSE),
  fit_imp3 = bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi(se) ~ age) +
    set_rescor(FALSE))
d <- nhanes
d$se <- se_col
for (nm in names(forms)) {
  w <- character()
  f <- withCallingHandlers(suppressMessages(frm(forms[[nm]], data = d)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  p <- f$opt$par
  se_fit <- suppressWarnings(sqrt(diag(vcov(f))))
  H0 <- stats::optimHess(p, f$obj$fn, f$obj$gr)
  Hr <- numDeriv::jacobian(function(x) as.numeric(f$obj$gr(x)), p)
  Hr <- (Hr + t(Hr)) / 2
  se0 <- sqrt(diag(solve(H0)))
  ser <- sqrt(diag(solve(Hr)))
  # vcov() names its rows by coefficient; the outer vector is in
  # template order
  tn <- gsub("[()]", "", names(par_template(f)$beta))
  se_fit <- se_fit[tn]
  k <- seq_along(se_fit)
  cat(sprintf("\n%s: code %d, warnings %d, logLik %.9f\n", nm,
              f$opt$convergence, length(w), as.numeric(logLik(f))))
  print(data.frame(par = names(se_fit), vcov = signif(se_fit, 6),
                   optimHess = signif(se0[k], 6),
                   richardson = signif(ser[k], 6),
                   vcov_rel = signif(se_fit / ser[k] - 1, 3),
                   optimHess_rel = signif(se0[k] / ser[k] - 1, 3)))
}

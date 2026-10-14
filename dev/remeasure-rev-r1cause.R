# Reviewer probe for R1: is 0.70.0's non-PD Hessian of fit_imp2 the
# optimHess step taken in par_units coordinates (autoscale_sdreport(),
# R/autoscale.R), rather than the objective?
#   Rscript dev/remeasure-rev-r1cause.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "\n")
data("nhanes", package = "mice")
bform2 <- bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
  set_rescor(FALSE)
f <- suppressWarnings(frm(bform2, data = nhanes))
p <- f$opt$par
u <- f$par_units %||% rep(1, length(p))
cat("par_units:", signif(u, 4), "\n")
ev <- function(H) signif(range(eigen((H + t(H)) / 2, symmetric = TRUE,
                                     only.values = TRUE)$values), 4)
for (nd in c(1e-3, 1e-4, 1e-5)) {
  Hq <- stats::optimHess(p / u, function(q) f$obj$fn(q * u),
                         function(q) f$obj$gr(q * u) * u,
                         control = list(ndeps = rep(nd, length(p))))
  H <- Hq / outer(u, u)
  se <- tryCatch(sqrt(diag(solve(H))), error = function(e) NA)
  cat(sprintf("q-space optimHess ndeps %g: eigen %s; se %s\n", nd,
              paste(ev(Hq), collapse = " to "),
              paste(signif(se, 4), collapse = " ")))
}
H0 <- stats::optimHess(p, f$obj$fn, f$obj$gr)
cat("raw-space optimHess ndeps 1e-3: eigen", ev(H0), "; se",
    signif(sqrt(diag(solve(H0))), 4), "\n")
# the inner problem: does the inner Newton reach its tolerance at the
# optimum with this build's inner.control?
cat("inner.control:\n"); str(f$obj$env$inner.control)

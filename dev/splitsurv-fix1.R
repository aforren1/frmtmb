f <- "R/curve.R"
x <- readLines(f)
rep1 <- function(old, new) {
  i <- which(x == old)
  if (length(i) != 1L) stop("not found once: ", old)
  x[i] <<- new
}
rep1("#' `gratia::confint(type = \"simultaneous\")` uses on an mgcv fit, and this",
     "#' `gratia::confint(type = \"simultaneous\")` uses on an mgcv fit, and")
rep1("#' package's simulation reproduces gratia's critical value inside its",
     "#' this simulation reproduces gratia's critical value inside its")
i <- grep("^#' Up to frmtmb 0.51.0 there was no such seam", x)
x <- c(x[1:(i - 1)],
  "#' Up to frmtmb 0.51.0 there was no such seam. frmtmb.spline, where",
  "#' this function lived until it moved into frmtmb, rebuilt `A` by unit",
  "#' perturbation, one `predict()` call per contributing coefficient, and",
  "#' read `V` out of `fit$cache$Vjoint`, which was an internal with no",
  "#' precedent to point at. One call replaced the reconstruction and the",
  "#' reach, and what the covariance check verifies has changed from \"the",
  "#' reconstruction reproduced core's number\" to \"the seam is being read",
  "#' correctly\".",
  x[(i + 8):length(x)])
rep1("#' is NOT the difference between two new subjects, which this package",
     "#' is NOT the difference between two new subjects, which this function")
rep1("                          allow_new_levels, tol, contrast)",
     "                       allow_new_levels, tol, contrast)")
rep1("  curve_assemble(parts, parts$eta, parts$se, parts$Sigma, level, simultaneous,",
     "  curve_assemble(parts, parts$eta, parts$se, parts$Sigma, level,")
rep1("              nsim, transform, seed, newdata,",
     "                 simultaneous, nsim, transform, seed, newdata,")
rep1("              what = if (is.null(contrast)) \"value\" else \"difference\")",
     "                 what = if (is.null(contrast)) \"value\" else \"difference\")")
rep1("curve_assemble <- function(parts, est, se, Sigma, level, simultaneous, nsim,",
     "curve_assemble <- function(parts, est, se, Sigma, level, simultaneous,")
rep1("                        transform, seed, newdata, what) {",
     "                           nsim, transform, seed, newdata, what) {")
rep1("  lk <- curve_predict_eta(fit, nd, parts$dpar, parts$resp, parts$re_formula,",
     "  lk <- curve_predict_eta(fit, nd, parts$dpar, parts$resp,")
rep1("                       parts$allow_new_levels)",
     "                          parts$re_formula, parts$allow_new_levels)")
i <- which(x == "`%||%` <- function(a, b) if (is.null(a)) b else a")
x <- x[-c(i, i + 1L)]
writeLines(x, f)

# Punch round 1: fd_hessian() against stats::optimHess() on the same
# objective, bitwise, so the fit-time Hessian sdreport() is handed is
# the one it would build. Models with and without random effects and an
# autoscaled one.
#   Rscript dev/nanse-p1-fdident.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-nanse-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r5",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb from", find.package("frmtmb"), "\n")
set.seed(1)
sleep <- lme4::sleepstudy
epi <- brms::epilepsy
dd <- data.frame(x = rnorm(200) * 1e-5, z = rnorm(200))
dd$y <- 1 + 3e4 * dd$x + 0.5 * dd$z + rnorm(200)
fits <- list(
  lm = frm(Reaction ~ Days, data = sleep),
  lmm = frm(Reaction ~ Days + (Days | Subject), data = sleep),
  glmm = frm(count ~ zAge + zBase * Trt + (1 | patient), data = epi,
             family = poisson()),
  autoscaled = suppressWarnings(frm(y ~ x + z, data = dd,
                                    control = frmtmb_control(autoscale = TRUE)))
)
for (nm in names(fits)) {
  f <- fits[[nm]]
  obj <- f$obj
  u <- f$par_units
  st <- ns$obj_state_save(obj)
  if (is.null(u) || all(u == 1)) {
    at <- ns$sdr_outer_point(obj)
    a <- stats::optimHess(at, obj$fn, obj$gr)
    ns$obj_state_restore(obj, st)
    b <- ns$fd_hessian(at, obj$gr)$H
  } else {
    q0 <- f$opt$par / u
    a <- stats::optimHess(q0, function(q) obj$fn(q * u),
                          function(q) obj$gr(q * u) * u)
    ns$obj_state_restore(obj, st)
    b <- ns$fd_hessian(q0, function(q) obj$gr(q * u) * u)$H
  }
  ns$obj_state_restore(obj, st)
  cat(sprintf("%-11s p=%d identical %s  max |diff| %.3g\n", nm, nrow(a),
              identical(a, b), max(abs(a - b))))
}

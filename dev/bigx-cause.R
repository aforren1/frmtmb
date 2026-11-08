# Lane bigx: the intermediate quantities of poisson y ~ (z + m) + (1 | g)
# on one data seed, per offset m: where the optimizer starts and stops,
# its code, the gradient at the stop, how far the stop is below the
# maximum, whether the tape itself is right at the maximum (the inner
# solve), the condition number of the outer Hessian at the maximum, the
# design's condition number, whether the default autoscale engages, and
# whether the Newton polish's Hessian (optimHess(), absolute steps) is
# finite at the stop.
#   BIGX_ARM=before|after Rscript dev/bigx-cause.R [seed]
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/bigx-common.R")
ns <- asNamespace("frmtmb")
a <- commandArgs(TRUE)
s <- if (length(a)) as.integer(a[1]) else 2L
for (m in c(1e4, 1e5, 1e6, 1e7, 1e8)) {
  d <- bigx_data("poisson", m, s)
  fc <- frm(bf(y ~ z + (1 | g)), family = poisson(), data = d)
  pc <- fc$opt$par
  pr <- pc
  pr[1] <- pc[1] - pc[2] * m
  f <- suppressWarnings(suppressMessages(
    frm(bf(y ~ w + (1 | g)), family = poisson(), data = d)))
  cat(sprintf("m %g seed %d\n", m, s))
  if (length(f$opt$par) != 3L) {
    cat("  slope dropped by the frame; parameters:",
        paste(names(f$opt$par), collapse = ","), "\n")
    next
  }
  st <- f$obj$par
  cat("  start (obj$par):", format(st, digits = 6), "\n")
  cat("  stop:", format(f$opt$par, digits = 8), " code", f$opt$convergence,
      f$opt$message, " iterations", f$opt$iterations, "\n")
  g <- as.numeric(f$obj$gr(f$opt$par))
  cat("  gradient at stop:", format(g, digits = 3), "\n")
  cat(sprintf("  objective: stop %.9f, raw tape at the centered maximum %.9f, centered fit %.9f; short by %.4g\n",
              f$opt$objective, f$obj$fn(pr), fc$opt$objective,
              f$opt$objective - fc$opt$objective))
  steps <- ns$outer_fd_steps(f$frame, names(pr))
  H <- ns$fd_hessian(pr, function(q) as.numeric(f$obj$gr(q)), eps = steps)$H
  D <- sqrt(abs(diag(H)))
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  cat("  outer Hessian at the maximum (fd, outer_fd_steps): unit-diagonal",
      "condition", format(max(abs(ev)) / min(abs(ev)), digits = 3),
      " smallest/largest", format(min(ev) / max(ev), digits = 3), "\n")
  cat("  design kappa", format(ns$predictor_design_kappa(f$frame), digits = 3),
      " default autoscale plan:",
      if (is.null(ns$autoscale_decide(f$frame, frmtmb_control(), FALSE)))
        "none" else "engaged",
      " column sd", format(sd(d$w), digits = 3), "\n")
  Ho <- tryCatch(optimHess(f$opt$par, f$obj$fn, f$obj$gr),
                 error = function(e) NULL)
  cat("  optimHess() at the stop: non-finite entries",
      if (is.null(Ho)) "error" else sum(!is.finite(Ho)), "of 9\n")
  cat("  optimizer basis on the objective:",
      !is.null(f$obj$env$.frmtmb_opt_basis), "\n")
}

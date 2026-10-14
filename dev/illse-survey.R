# Lane illse: which fits of the design families get which SE verdict,
# with the condition number of the outer Hessian scaled to unit
# diagonal (the AD Hessian without random effects, the finite-difference
# one with them). Seeds 1 to 3.
#   ILLSE_ARM=base|lane Rscript dev/illse-survey.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-designs.R")
tag <- function(w) {
  if (!length(w)) return("-")
  k <- c(illcond = "ill-conditioned", short = "stopped short",
         lost = "Standard errors are not available",
         conv = "onvergence", grad = "gradient", flat = "flat",
         boundary = "boundary", sep = "separat")
  paste(names(k)[vapply(k, function(s) any(grepl(s, w)), NA)],
        collapse = ",")
}
kap <- function(f) {
  H <- tryCatch(if (length(f$obj$env$random)) {
    frmtmb:::fit_outer_hessian(f)$H
  } else f$obj$he(f$opt$par), error = function(e) NULL)
  if (is.null(H) || !all(is.finite(H))) return(NA_real_)
  D <- sqrt(abs(diag(H)))
  ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
  max(abs(ev)) / min(ev)
}
for (mk in illse_designs()) for (seed in 1:3) {
  ds <- mk(seed)
  t0 <- proc.time()[[3]]
  r <- cap(frm(ds$formula, data = ds$data, family = ds$family))
  el <- proc.time()[[3]] - t0
  if (inherits(r$v, "illse_err")) {
    cat(sprintf("%-34s s%d ERROR %s\n", ds$lab, seed, substr(r$v, 1, 70)))
    next
  }
  f <- r$v
  se <- tryCatch(suppressWarnings(sqrt(diag(vcov(f)))),
                 error = function(e) NA)
  cat(sprintf("%-34s s%d code %d | lost %-14s | warn %-26s | kappa %.1e | finite SE %d/%d | %.1fs\n",
              ds$lab, seed, f$opt$convergence, lost_of(f), tag(r$w),
              kap(f), sum(is.finite(se)), length(se), el))
}

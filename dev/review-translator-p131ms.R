# Reviewer scratch (review of lane translator, 2026-11-09): p131
# (cs_pred() x mixture). The record's L4 route left frm()'s maximum 0.814
# nats BELOW brms's best of 6 optimizing() starts. With the lane's L1
# identity the two log densities are the same function, so either brms's
# point is a better point of frm()'s own objective (frm() stopped short)
# or the gap comes from elsewhere. frm() is run with an optimizer that
# tries 20 seeded random starts around frm()'s own start and keeps the
# best. Records nothing.
#
#   Rscript dev/review-translator-p131ms.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
e <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = e)
}
spec <- e$eq_specs[["cs_pred() x mixture"]]()
f0 <- suppressWarnings(e$eq_frm_fit(spec))
cat(sprintf("frm() default: logLik %.6f\n", as.numeric(logLik(f0))))
print(f0$estimates)
ms <- function(par, fn, gr, lower, upper, control) {
  best <- NULL
  vals <- numeric()
  for (s in 0:20) {
    set.seed(1000 + s)
    p0 <- if (s == 0) par else par + rnorm(length(par), 0, 1)
    p0 <- pmin(pmax(p0, lower), upper)
    r <- tryCatch(stats::nlminb(p0, fn, gr, lower = lower, upper = upper,
                                control = list(iter.max = 2000,
                                               eval.max = 4000)),
                  error = function(err) NULL)
    if (is.null(r) || !is.finite(r$objective)) {
      vals <- c(vals, NA)
      next
    }
    vals <- c(vals, -r$objective)
    if (is.null(best) || r$objective < best$objective) best <- r
  }
  cat("starts (minus objective):", paste(round(vals, 4), collapse = " "),
      "\n")
  list(par = best$par, objective = best$objective,
       convergence = best$convergence, message = best$message)
}
f1 <- suppressWarnings(e$eq_frm_fit(spec, ctrl = list(optimizer = ms,
                                                      restarts = 0,
                                                      check_se = "ignore")))
cat(sprintf("frm() best of 21 starts: logLik %.6f\n",
            as.numeric(logLik(f1))))
print(f1$estimates)

# Review of lane fitedge: the theta layout of each fresh design, and one fit.
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs("lane")
source(file.path(FE_WT, "dev/review-fitedge-designs.R"))
for (ds in rv_designs) {
  d <- rv_mk(ds, 1, 60, 1)
  sp <- rv_spell[[ds]]
  t0 <- proc.time()[[3]]
  r <- fe_fit(sp$nl, data = d, family = sp$fam)
  if (inherits(r$fit, "error")) {
    cat(ds, "ERROR", conditionMessage(r$fit), "\n")
    next
  }
  th <- r$fit$opt$par[names(r$fit$opt$par) == "theta"]
  cat(ds, "ll", format(as.numeric(logLik(r$fit)), digits = 10), "conv",
      r$fit$opt$convergence, "ntheta", length(th), "theta",
      paste(round(th, 3), collapse = " "), "sec",
      round(proc.time()[[3]] - t0, 2), "warn", length(r$warn), "\n")
  print(r$fit$frame$re_blocks |> lapply(function(b) b[c("covstruct")]) |>
          unlist())
  if (length(r$warn)) cat("  W:", substr(r$warn, 1, 120), sep = "\n  ")
  cat("\n")
}

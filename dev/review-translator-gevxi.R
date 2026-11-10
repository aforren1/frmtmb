# Reviewer scratch (review of lane translator, 2026-11-09): why
# dev/translator-gev-check.R warns "longer object length is not a multiple
# of shorter object length": how many distinct values do xi_k and sigma_k
# take over the rows at frm()'s estimate for p021? Records nothing.
.libPaths(c("C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
e <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = e)
}
spec <- e$eq_specs[["gen_extreme_value x mixture"]]()
fit <- e$eq_frm_fit(spec)
dp <- frmtmb:::eval_dpars(fit)[[1L]]
for (nm in c("sigma1", "xi1", "sigma2", "xi2")) {
  v <- as.numeric(dp[[nm]])
  u <- unique(v)
  cat(sprintf("%-6s length %d, distinct %d, range %.17g .. %.17g\n", nm,
              length(v), length(u), min(v), max(v)))
}

# Review of lane sampfront: prior_samples(), model_weights(), pp_average()
# reach brms's methods for a brmsfit when frmtmb.sample is attached after
# brms, and reach frmtmb.sample's methods when attached before.
#   Rscript dev/review-sampfront-owner.R [after|before]
ord <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-sampfront-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  if (ord == "after") { library(brms); library(frmtmb.sample) } else {
    library(frmtmb.sample); library(brms) }
})
cat("order: frmtmb.sample", ord, "brms\n")
for (f in c("prior_samples", "model_weights", "pp_average", "prior_draws")) {
  cat(f, "found in:", environmentName(environment(get(f))), "\n")
}
fake <- structure(list(), class = "brmsfit")
for (f in c("prior_samples", "model_weights", "pp_average")) {
  r <- tryCatch(withCallingHandlers(get(f)(fake),
                  warning = function(w) {
                    cat("  warning:", conditionMessage(w), "\n")
                    invokeRestart("muffleWarning")
                  }), error = function(e) e)
  cat(f, "(brmsfit):", if (inherits(r, "error")) paste("ERROR in",
      deparse(conditionCall(r))[1], ":", substr(conditionMessage(r), 1, 90))
      else "returned", "\n")
}

## Reviewer: fit the review models in frmtmb (lane build) and sample them
## with frmtmb.sample, sampler seed 3. Saves fits and draws.
##   Rscript dev/brmsnames-rev-frm.R
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
source("dev/brmsnames-rev-data.R")
ms <- rev_models("frmtmb")
av <- commandArgs(trailingOnly = TRUE)
if (length(av)) ms <- ms[av]
for (nm in names(ms)) {
  out <- sprintf("dev/stan-cache/brmsnames-rev-frm-%s.rds", nm)
  fit <- tryCatch(q(frm(ms[[nm]]$f, family = gaussian(), data = rev_dd,
                        start = if (nm == "C5") list(beta = c(2, 0.7)))),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    cat(nm, "frm ERROR", conditionMessage(fit), "\n"); next
  }
  ds <- tryCatch(q(frm_sample(fit, chains = 2, iter = 200, refresh = 0,
                              seed = 3)), error = function(e) e)
  if (inherits(ds, "error")) cat(nm, "sample ERROR", conditionMessage(ds), "\n")
  saveRDS(list(fit = fit, ds = ds), out)
  cat(nm, "variables(fit):", variables(fit), "\n")
  if (!inherits(ds, "error")) cat(nm, "draws cols:", head(colnames(ds$draws), 40), "\n")
}
cat("DONE\n")

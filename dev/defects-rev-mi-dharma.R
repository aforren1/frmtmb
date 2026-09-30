# Reviewer of lane defects, recheck: dharma_residuals() and osa residuals
# on the univariate mi() fits, per arm.  Rscript ... lane|base
arm <- commandArgs(trailingOnly = TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-defects-lib", libs)
.libPaths(libs)
suppressMessages(library(frmtmb))
source("dev/defects-rev-mi-data.R")
d <- mi_data()
m <- function(e) tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
for (nm in c("uni", "sdy")) {
  f <- frm(eval(mi_models[[nm]]), data = d)
  set.seed(3)
  dh <- m(suppressWarnings(dharma_residuals(f)))
  cat(nm, "dharma:", if (is.character(dh)) dh else
      c(n_obs = length(dh$observedResponse), nsim_rows = nrow(dh$simulatedResponse),
        n_fitted = length(dh$fittedPredictedResponse),
        placeholder_zero_in_obs = sum(dh$observedResponse == 0)), "\n")
  cat(nm, "osa:", substr(m(residuals(f, type = "osa")), 1, 120), "\n")
}

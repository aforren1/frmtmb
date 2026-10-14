.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
dd <- fs_data_dirichlet(seed = 11)
f1 <- frm(bf(Y ~ x), family = dirichlet(), data = dd)
for (ty in c("osa", "deviance", "quantile")) {
  cat("--- residuals type", ty, "\n")
  cat(tryCatch({r <- residuals(f1, type = ty); paste(class(r), dim(r))},
               error = function(e) conditionMessage(e)), "\n")
}
cat("--- dharma_residuals\n")
r <- tryCatch(dharma_residuals(f1, n = 50, seed = 1), error = function(e) conditionMessage(e))
if (is.character(r)) cat(r, "\n") else { cat(class(r), "\n"); str(r$scaledResiduals) }
cat("--- diagnose\n")
r <- tryCatch(capture.output(diagnose(f1)), error = function(e) conditionMessage(e))
cat(head(r, 30), sep = "\n")

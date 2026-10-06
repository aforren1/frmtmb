# The after-table of the user's decision of 2026-10-06 on predict() at
# crossing cs() rows, on the release review's fit and grid
# (dev/relrev-cs.R, dev/relrev-log/cs-fit.rds): per row, the NA
# replicates of 1000, whether fitted() is NaN, and predict()'s
# proportions.
#   Rscript dev/rel068-m3-table.R > dev/rel068-log/m3-table.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
S <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev-log/cs-fit.rds")
set.seed(11)
raw <- predict(S$fit, newdata = S$nd, summary = FALSE)
set.seed(11)
p <- withCallingHandlers(predict(S$fit, newdata = S$nd), warning = function(w) {
  cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
})
fe <- fitted(S$fit, newdata = S$nd)[, "Estimate", ]
tab <- data.frame(row = seq_len(nrow(S$nd)), x = round(S$nd$x, 3),
                  na_reps = colSums(is.na(raw)),
                  fitted_nan = rowSums(is.nan(fe)) > 0, round(p, 4),
                  check.names = FALSE)
print(tab, row.names = FALSE)

# Review probe: what a read of wiener_gng()'s CDF slot by `$` gave on
# 0.83.0 (rellib-r22) and gives on the lane's build.
#   Rscript dev/reviews/cdfname-review-gng-base.R [base]
args <- commandArgs(trailingOnly = TRUE)
base <- length(args) && args[[1L]] == "base"
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-cdfname-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.eam"), "\n")
slot <- if (base) "lcdf" else "cdf"
f <- wiener_gng(deadline = 1.5)
r <- tryCatch({
  v <- f[[slot]]
  w <- do.call(`$`, list(f, slot))
  sprintf("[[ ]] %s, $ %s", if (is.null(v)) "NULL" else "set",
          if (is.null(w)) "NULL" else "set")
}, error = function(e) paste("ERROR", conditionMessage(e)))
cat("G wiener_gng() slot", slot, ":", r, "\n")

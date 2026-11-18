# rel084-kt copy of dev/reviews/slopetrig-review-probe.R: lib = rellib-r23 then user lib (lane lib and rellib-r22 removed). No other change.
# Review of lane slopetrig: what ratio autoscale_plan() reads, and
# whether the default engages, on designs the lane did not calibrate on.
# One data set per design (seed 1); no fit beyond the default.
# Usage: Rscript slopetrig-review-probe.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r23",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
source("dev/reviews/slopetrig-review-designs.R")
for (nm in names(designs)) {
  ds <- designs[[nm]]
  d <- ds$data(1)
  r <- tryCatch({
    f <- suppressWarnings(suppressMessages(
      frm(ds$form, data = d, family = ds$family)))
    pl <- ns$autoscale_plan(f$frame)
    zs <- unlist(lapply(pl, function(p) lapply(p$z, function(z)
      sprintf("%s ratio=%.3g spread=%.3g", z$label, z$ratio, z$spread))))
    paste0("engaged=", isTRUE(f$autoscaled), " | ",
           paste(zs, collapse = " ; "))
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-14s %s\n", nm, r))
}

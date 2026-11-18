# Review probe: does frmtmb.eam's ndt seam drop the `cdf` element of a
# family that has none, so that `$cdf` on the finalized family errors?
#   Rscript dev/reviews/cdfname-review-eam.R
.libPaths(c("C:/Users/adf44/source/r/wt-cdfname-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.eam"), "\n")
say <- function(label, expr) {
  r <- tryCatch({
    v <- force(expr)
    if (is.null(v)) "value NULL" else paste("value of class", class(v)[1L])
  }, error = function(e) sprintf("ERROR [%s] %s", class(e)[1L],
                                 substr(conditionMessage(e), 1, 120)))
  cat(sprintf("E %-40s %s\n", label, r))
}
for (nm in c("wiener_gng", "gddm", "lba")) {
  f <- switch(nm, wiener_gng = wiener_gng(deadline = 1.5), gddm = gddm(),
              lba = lba(2))
  cat("E", nm, "has cdf element:", "cdf" %in% names(f), "\n")
  fin <- tryCatch(frmtmb.eam:::ddm_ndt_install(f, ub = 0.2, floors = NULL,
                                               what = "rt"),
                  error = function(e) {
                    cat("E install error:", conditionMessage(e), "\n")
                    NULL
                  })
  if (!is.null(fin)) {
    cat("E", nm, "after ndt install has cdf element:",
        "cdf" %in% names(fin), "\n")
    say(paste0(nm, " finalized $cdf"), fin$cdf)
    say(paste0(nm, " finalized [[\"cdf\"]]"), fin[["cdf"]])
  }
}
say("wiener()$lcdf", wiener()$lcdf)

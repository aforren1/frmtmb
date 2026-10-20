# Lane sampleprior: which families frmtmb and brms each name, and the
# dpars of each. Usage: Rscript dev/sampleprior-fams.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[nzchar(a)], "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(brms)})
cat("lib:", find.package("frmtmb"), "\n")
reg <- sort(names(frmtmb:::family_registry))
cat("frmtmb registry:", paste(reg, collapse = ", "), "\n")
for (nm in reg) {
  f <- tryCatch(frmtmb::brmsfamily(nm), error = function(e) NULL)
  if (is.null(f)) f <- tryCatch(frmtmb:::as_frmtmb_family(nm),
                                error = function(e) NULL)
  if (is.null(f)) { cat(nm, ": (no constructor)\n"); next }
  lk <- vapply(f$links %||% list(), function(l) l$name %||% "?", "")
  cat(sprintf("%-28s %s\n", nm, paste0(names(lk), "[", lk, "]",
                                         collapse = " ")))
}

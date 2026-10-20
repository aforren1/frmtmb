# Lane prilabel: how many standard deviations each covariance structure
# carries at d = 3, which decides whether a per-coefficient sd prior is
# one parameter of it.
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
reg <- frmtmb:::covstruct_registry
for (nm in names(reg)) {
  f <- reg[[nm]]$sd_idx
  v <- if (is.function(f)) tryCatch(paste(f(3L), collapse = ","),
                                    error = function(e) "ERR") else "none"
  np <- tryCatch(reg[[nm]]$npar(3L), error = function(e) NA)
  cat(sprintf("%-10s sd_idx(3) = %-10s npar(3) = %s\n", nm, v, np))
}

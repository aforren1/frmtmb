# pkgdown::check_pkgdown() for all nine packages on the release tree.
#   Rscript dev/rel075-pkgdown.R > dev/rel075-log/pkgdown.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-rel075"
pk <- c("frmtmb", paste0("frmtmb.", c("coupling", "eam", "latent", "learn",
                                      "ode", "sample", "spatial", "survival")))
for (p in pk) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  r <- tryCatch({ pkgdown::check_pkgdown(path); "OK" },
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("PKGDOWN", p, r, "\n")
}

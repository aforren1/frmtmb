# pkgdown::check_pkgdown() for all eight packages on the release tree.
#   Rscript dev/rel069-pkgdown.R > dev/rel069-log/pkgdown.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
pk <- c("frmtmb", paste0("frmtmb.", c("coupling", "eam", "latent", "learn",
                                      "ode", "sample", "spline")))
for (p in pk) {
  path <- if (p == "frmtmb") root else file.path(root, "extensions", p)
  r <- tryCatch({ pkgdown::check_pkgdown(path); "OK" },
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("PKGDOWN", p, r, "\n")
}

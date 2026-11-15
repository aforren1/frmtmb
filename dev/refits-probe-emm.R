.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (p in c("emmeans", "marginaleffects", "brms", "insight", "data.table",
            "rstan", "mgcv")) {
  cat(p, as.character(packageVersion(p)), "\n")
}
f <- deparse(emmeans::ref_grid)
i <- grep("recover_data|try[(]|is.character[(]data|emm_basis", f)
for (k in i) cat(k, ":", f[(k):(k + 2)], sep = "\n")

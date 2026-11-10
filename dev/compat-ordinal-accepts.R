# Which addition terms the ordinal families declare, and which the
# custom family path declares, on the library named by the first
# argument (default: the base build).
#   Rscript dev/compat-ordinal-accepts.R [lib;lib...]
args <- commandArgs(trailingOnly = TRUE)
lib <- if (length(args) >= 1) args[1] else
  "C:/Users/adf44/source/r/rellib-r19"
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
ns <- asNamespace("frmtmb")
acc <- ns$compat_core_family_accepts()
ft <- ns$frmtmb_compat_features_tbl()
cat("aterm features:", paste(ft$name[ft$kind == "aterm"], collapse = " "),
    "\n")
for (f in c(ns$frmtmb_compat_groups_lst$ordinal, "gaussian", "poisson")) {
  a <- acc[[f]]
  cat(sprintf("%-18s %s\n", f,
              if (is.null(a)) "<NULL: no allow-list>" else
                paste(a, collapse = " ")))
}

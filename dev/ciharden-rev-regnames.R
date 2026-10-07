# Reviewer: which call-shaped registry names are functions of a
# base-priority package that is not loaded at startup?
.libPaths(c("C:/Users/adf44/source/r/wt-ciharden-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb.sample))
ft <- frmtmb::frm_compat_features()
nm <- sub("[(][)]$", "", grep("[(][)]$", ft$name, value = TRUE))
ref <- frmtmb::frm_compat("frm_sample", status = "refused")
cat("call-shaped names:", length(nm), "\n")
cat("refused for frm_sample:", nrow(ref), "rows\n")
pk <- c("splines", "stats4", "grid", "tools", "parallel", "compiler",
        "tcltk")
for (p in pk) {
  ns <- asNamespace(p)
  hit <- nm[vapply(nm, function(n) {
    exists(n, ns, inherits = FALSE) && is.function(get(n, ns))
  }, NA)]
  cat(p, ":", if (length(hit)) paste(hit, collapse = " ") else "-", "\n")
}

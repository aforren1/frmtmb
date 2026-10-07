# base_r_function() decides an unloaded base-priority package by
# exists() alone. Does any call-shaped registry name exist as an object
# (function or not) in one of the packages a session may not load?
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ft <- frmtmb::frm_compat_features()
nm <- sub("[(][)]$", "", grep("[(][)]$", ft$name, value = TRUE))
cat("call-shaped names:", paste(nm, collapse = " "), "\n")
for (p in c("splines", "stats4", "grid", "tools", "parallel", "compiler",
            "tcltk")) {
  e <- new.env(parent = emptyenv())
  lazyLoad(file.path(find.package(p, lib.loc = .Library), "R", p),
           envir = e)
  hit <- nm[vapply(nm, exists, NA, envir = e, inherits = FALSE)]
  cat(p, ":", if (length(hit)) paste(hit, collapse = " ") else "-", "\n")
}

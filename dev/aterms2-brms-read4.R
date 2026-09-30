.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
fns <- ls(ns, all.names = TRUE)
for (nm in c(".frame_index", "vars_keep_na.brmsterms",
             "vars_keep_na.mvbrmsterms", "get_ad_vars")) {
  cat("\n######## ", nm, "\n")
  if (exists(nm, envir = ns)) print(get(nm, envir = ns))
}
cat("\n######## families accepting rate / subset / index / cat / thres\n")
famnames <- sub("^\\.family_", "", grep("^\\.family_", fns, value = TRUE))
for (f in famnames) {
  ad <- tryCatch(get(paste0(".family_", f), envir = ns)()$ad,
                 error = function(e) NULL)
  cat(sprintf("%-28s %s\n", f,
              paste(intersect(ad, c("rate", "subset", "index", "cat",
                                    "thres", "mi", "trunc", "cens",
                                    "weights")),
                    collapse = " ")))
}

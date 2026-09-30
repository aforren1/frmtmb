.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
fns <- ls(ns)
hit <- function(pat) {
  out <- character()
  for (f in fns) {
    obj <- get(f, envir = ns)
    if (!is.function(obj)) next
    txt <- paste(deparse(obj), collapse = "\n")
    if (grepl(pat, txt)) out <- c(out, f)
  }
  out
}
cat("index frame:", hit("frame\\$index|\\$index <-|index\\]\\] <-"), "\n")
cat("cross:", hit("check_cross_formula_indexing"), "\n")
show <- function(nm) {
  cat("\n######## ", nm, "\n")
  print(get(nm, envir = ns))
}
for (nm in c("check_cross_formula_indexing", "frame_index", "vars_keep_na",
             "subset_data", "frame_sp")) {
  if (exists(nm, envir = ns)) show(nm) else cat("\n## missing:", nm, "\n")
}
cat("\n######## families accepting rate / subset / index / cat / thres\n")
famnames <- sub("^\\.family_", "", grep("^\\.family_", fns, value = TRUE))
for (f in famnames) {
  ad <- tryCatch(get(paste0(".family_", f), envir = ns)()$ad,
                 error = function(e) NULL)
  cat(sprintf("%-28s %s\n", f,
              paste(intersect(ad, c("rate", "subset", "index", "cat",
                                    "thres", "mi", "trunc", "cens")),
                    collapse = " ")))
}

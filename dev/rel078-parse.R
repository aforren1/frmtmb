# Every R file of the merged tree parses, and core R/ has no top-level
# name defined twice (a merge can keep both copies of a function).
#   Rscript dev/rel078-parse.R
setwd("C:/Users/adf44/source/r/frmtmb-wt-rel078")
fs <- c(list.files("R", "[.]R$", full.names = TRUE),
        list.files("extensions", "[.]R$", recursive = TRUE,
                   full.names = TRUE),
        list.files("tests", "[.]R$", recursive = TRUE, full.names = TRUE),
        list.files(".github", "[.]R$", recursive = TRUE, full.names = TRUE))
bad <- 0
for (f in fs) {
  r <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(r)) { bad <- bad + 1; cat("PARSE ERROR", f, r, "\n") }
}
cat("parsed", length(fs), "files, errors", bad, "\n")
for (d in c("R", file.path("extensions", list.files("extensions"), "R"))) {
  if (!dir.exists(d)) next
  defs <- unlist(lapply(list.files(d, "[.]R$", full.names = TRUE),
                        function(f) {
    e <- parse(f)
    vapply(e, function(x) {
      if (is.call(x) && as.character(x[[1]])[1] %in% c("<-", "=") &&
            is.name(x[[2]])) as.character(x[[2]]) else NA_character_
    }, "")
  }))
  defs <- defs[!is.na(defs)]
  cat(d, ": top-level names", length(defs), ", defined twice:",
      paste(unique(defs[duplicated(defs)]), collapse = ", "), "\n")
}

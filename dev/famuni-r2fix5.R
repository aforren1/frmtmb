p <- "R/families-extreme.R"
x <- paste(readLines(p), collapse = "\n")
old <- "  if (length(out) > 1L) {
    out <- lapply(out, function(tm) {
      tm$blocked <- c(tm$blocked, \"a second gen_extreme_value() response\")
      tm
    })
  }"
new <- "  # the REPORTed name is one per model, so two charted responses would
  # share it; a mixture's components are not charted anyway
  plain <- vapply(out, function(tm) !\"a mixture()\" %in% tm$blocked, NA)
  if (sum(plain) > 1L) {
    out[plain] <- lapply(out[plain], function(tm) {
      tm$blocked <- c(tm$blocked, \"a second gen_extreme_value() response\")
      tm
    })
  }"
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)

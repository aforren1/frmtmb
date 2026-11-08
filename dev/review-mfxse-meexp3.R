# Review of lane mfxse: marginaleffects 1.0.0 settings_get(), and where
# the jacobian_function setting is documented.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
ns <- asNamespace("marginaleffects")
print(get("settings_get", ns))
db <- tools::Rd_db("marginaleffects")
for (n in names(db)) {
  t <- paste(capture.output(tools::Rd2txt(db[[n]], options = list(underline_titles = FALSE))), collapse = "\n")
  if (grepl("jacobian_function", t)) {
    cat("=== in", n, "\n")
    l <- strsplit(t, "\n")[[1]]; i <- grep("jacobian_function", l)
    for (k in i) cat(l[max(1, k - 3):min(length(l), k + 6)], sep = "\n")
  }
  if (grepl("hypothesis_is_formula", t)) cat("hypothesis_is_formula documented in", n, "\n")
}

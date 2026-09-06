ns <- asNamespace("brms")
fns <- grep("^[.]family_", ls(ns, all.names = TRUE), value = TRUE)
rows <- list()
for (f in fns) {
  fam <- sub("^[.]family_", "", f)
  info <- tryCatch(get(f, envir = ns)(), error = function(e) NULL)
  if (is.null(info)) next
  lk <- info$links
  rows[[fam]] <- lk
  dp <- info$dpars
  cat(fam, "|", paste(lk, collapse = ","), "|dpars:", paste(dp, collapse = ","), "\n")
}
cat("=====dpar links=====\n")
dfns <- grep("^[.]dpar_", ls(ns, all.names = TRUE), value = TRUE)
for (f in dfns) {
  info <- tryCatch(get(f, envir = ns)(), error = function(e) NULL)
  if (is.null(info)) next
  cat(sub("^[.]dpar_", "", f), "|", paste(info$links, collapse = ","), "\n")
}
saveRDS(rows, "dev/lk-brms-famlinks.rds")

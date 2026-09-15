# After an interrupted session: every R file this lane touched parses,
# and the installed packages are newer than their newest source edit.
root <- "C:/Users/adf44/source/r/frmtmb-wt-samplegen"
fs <- system2("git", c("-C", root, "status", "--porcelain"), stdout = TRUE)
fs <- trimws(substring(fs, 4))
fs <- unlist(lapply(file.path(root, fs), function(p)
  if (dir.exists(p)) list.files(p, "[.]R$", recursive = TRUE, full.names = TRUE)
  else p))
fs <- fs[grepl("[.]R$", fs)]
bad <- 0L
for (f in fs) {
  ok <- tryCatch({ parse(f); TRUE }, error = function(e) {
    cat("PARSE FAIL", f, conditionMessage(e), "\n"); FALSE })
  if (!ok) bad <- bad + 1L
}
cat("R files parsed:", length(fs), " failures:", bad, "\n")
lib <- "C:/Users/adf44/source/r/samplegen-lib"
for (p in c("frmtmb", "frmtmb.sample")) {
  src <- if (p == "frmtmb") root else file.path(root, "extensions/frmtmb.sample")
  srcs <- c(list.files(file.path(src, "R"), full.names = TRUE),
            file.path(src, "NAMESPACE"))
  inst <- file.info(file.path(lib, p, "Meta", "package.rds"))$mtime
  newest <- max(file.info(srcs)$mtime)
  cat(p, "installed", format(inst), "newest source", format(newest),
      if (inst > newest) "OK" else "STALE", "\n")
}

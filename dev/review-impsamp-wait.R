# Wait until the named mutants have a RESULT line, up to about nine minutes.
ms <- commandArgs(TRUE)
d <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/mut2"
for (i in 1:54) {
  done <- vapply(ms, function(m) {
    f <- file.path(d, paste0("test-", m, ".log"))
    file.exists(f) && any(grepl("^RESULT|^Error", readLines(f, warn = FALSE)))
  }, TRUE)
  if (all(done)) break
  Sys.sleep(10)
}
print(done)

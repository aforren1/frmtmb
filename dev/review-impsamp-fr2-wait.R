# Wait up to 150 s for a line in the run's output, then print the stages reached.
a <- commandArgs(TRUE)
f <- "C:/Users/adf44/AppData/Local/Temp/2/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/tasks/b0jq02bo3.output"
for (i in 1:15) {
  l <- if (file.exists(f)) readLines(f, warn = FALSE) else character(0)
  if (a[[1]] %in% l) break
  Sys.sleep(10)
}
cat("stages:", l, "\n")

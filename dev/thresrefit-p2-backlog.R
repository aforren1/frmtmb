## Splice this lane's backlog entry into dev/test-backlog.md just above
## "## Reference", keeping the file's own line endings. Idempotent: an
## existing wt-thresrefit entry is cut out first.
doc <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit/dev/test-backlog.md"
add <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
              "c--Users-adf44-source-r-frmtmb/",
              "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/",
              "thresrefit-backlog.md")
old <- sub("\r$", "", readLines(doc, warn = FALSE))
i <- grep("^## Filed by wt-thresrefit", old)
j <- grep("^## Reference", old)
stopifnot(length(j) == 1L)
if (length(i)) {
  stopifnot(length(i) == 1L, i < j)
  old <- c(old[seq_len(i - 1L)], old[j:length(old)])
  j <- grep("^## Reference", old)
}
new <- sub("\r$", "", readLines(add, warn = FALSE))
writeLines(c(old[seq_len(j - 1L)], new, old[j:length(old)]), doc)
cat("spliced", length(new), "lines above line", j, "\n")

## Splice the punch-round section into dev/thresrefit-findings.md just
## above "## NEWS entry", keeping the file's own line endings. Idempotent:
## an existing punch-round section is cut out first.
doc <- "C:/Users/adf44/source/r/frmtmb-wt-thresrefit/dev/thresrefit-findings.md"
add <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
              "c--Users-adf44-source-r-frmtmb/",
              "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/",
              "thresrefit-p1-findings.md")
old <- sub("\r$", "", readLines(doc, warn = FALSE))
i <- grep("^## Punch round 1", old)
j <- grep("^## NEWS entry", old)
stopifnot(length(j) == 1L)
if (length(i)) {
  stopifnot(length(i) == 1L, i < j)
  old <- c(old[seq_len(i - 1L)], old[j:length(old)])
  j <- grep("^## NEWS entry", old)
}
new <- sub("\r$", "", readLines(add, warn = FALSE))
writeLines(c(old[seq_len(j - 1L)], new, old[j:length(old)]), doc)
cat("spliced", length(new), "lines above line", j, "\n")

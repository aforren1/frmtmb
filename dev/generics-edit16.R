root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
md <- file.path(root, "dev", "generics-findings.md")
blk <- readLines(file.path(root, "dev", "generics-out", "blocks",
                           "suite-final.txt"), warn = FALSE)
txt <- paste(readLines(md, warn = FALSE), collapse = "\n")
old <- paste0(
  "The core suite ran against the build in place before two\n",
  "documentation-only edits (shortening the markdown table cells in\n",
  "`R/scales.R` to 80 columns, and the comment blocks in\n",
  "`R/generic-owners.R`) and before the hook-idempotence fix. The two new\n",
  "test files were re-run against the FINAL build, 23 of 23 and 25 of 25.")
new <- paste0(
  "That run was against the build in place before two\n",
  "documentation-only edits and before the hook-idempotence fix, so the\n",
  "whole suite was run again, one process per file, against the FINAL\n",
  "build:\n\n",
  paste(blk, collapse = "\n"), "\n\n",
  "The two runs agree to the assertion: 1454 blocks and 10176\n",
  "assertions both times, 0 failures and 0 errors both times. The two\n",
  "new test files were also re-run against the final build on their own,\n",
  "23 of 23 and 25 of 25.")
if (!grepl(old, txt, fixed = TRUE)) stop("no match")
writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                    fixed = TRUE)[[1]], md, useBytes = TRUE)
cat("DONE\n")

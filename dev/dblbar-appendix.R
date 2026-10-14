# Lane dblbar: every test file's RESULT line, by tier and package, from
# the tier logs in dev/dblbar-tests/, for the findings' appendix.
# Output: dev/dblbar-tests/appendix.md
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
tiers <- c(final_all = "final-all", gated = "lane-gated",
           drmtmb = "lane-drm", last_edits = "lane-last")
out <- character(0)
for (k in names(tiers)) {
  lg <- readLines(file.path(wt, "dev/dblbar-tests",
                            paste0(tiers[[k]], ".log")))
  out <- c(out, "", paste0("### ", tiers[[k]]), "", "```")
  pkg <- ""
  for (l in lg) {
    m <- regmatches(l, regexec("^== (\\S+) [0-9]+ files ==$", l))[[1]]
    if (length(m)) pkg <- m[2]
    if (startsWith(l, "RESULT ")) {
      out <- c(out, paste(pkg, sub("^RESULT ", "", l)))
    }
  }
  out <- c(out, "```")
}
writeLines(out, file.path(wt, "dev/dblbar-tests/appendix.md"))
cat(length(grep("^frmtmb", out)), "result lines\n")

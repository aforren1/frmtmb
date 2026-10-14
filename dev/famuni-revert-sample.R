# Rebuild frmtmb.sample's R/sample.R as it was BEFORE fix round 1, by
# running this round's edit scripts backwards on a copy, for the "seen to
# fail first" arm of the new sample tests. Each script edits with
# r(old, new); here r() swaps the two. Writes the copy to the path given.
# Usage: Rscript dev/famuni-revert-sample.R <out.R>
setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
out <- commandArgs(TRUE)[1]
src <- "extensions/frmtmb.sample/R/sample.R"
x <- paste(readLines(src), collapse = "\n")
undo <- function(script) {
  e <- new.env()
  # each script's own r() and fixed replacements, read and replayed in
  # reverse on x
  code <- parse(script)
  pairs <- list()
  e$r <- function(old, new, all = FALSE) {
    pairs[[length(pairs) + 1L]] <<- list(old = old, new = new, all = all)
  }
  for (ex in code) {
    txt <- paste(deparse(ex), collapse = "\n")
    if (grepl("^r[(]", txt)) eval(ex, e)
  }
  if (grepl("fix5e", script)) {
    # this one used sub() directly with old/new variables
    eval(code[[grep("^old <-", vapply(code, function(z) {
      paste(deparse(z), collapse = "")
    }, ""))[1]]], e)
    eval(code[[grep("^new <-", vapply(code, function(z) {
      paste(deparse(z), collapse = "")
    }, ""))[1]]], e)
    pairs <- list(list(old = e$old, new = e$new, all = FALSE))
  }
  for (pr in rev(pairs)) {
    if (!grepl(pr$new, x, fixed = TRUE)) {
      cat("  skip (not in this file):", substr(pr$new, 1, 50), "\n")
      next
    }
    x <<- if (isTRUE(pr$all)) gsub(pr$new, pr$old, x, fixed = TRUE) else
      sub(pr$new, pr$old, x, fixed = TRUE)
  }
}
for (s in c("dev/famuni-fix5g.R", "dev/famuni-fix5f.R", "dev/famuni-fix5e.R",
            "dev/famuni-fix5d.R", "dev/famuni-fix5c.R")) {
  cat(s, "\n")
  undo(s)
}
writeLines(strsplit(x, "\n")[[1]], out)
cat("wrote", out, "\n")

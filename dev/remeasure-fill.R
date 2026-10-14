# Fill the generated blocks of dev/remeasure-findings.md from the
# reducers' output, so no count in the record is typed by hand.
#
#   Rscript dev/remeasure-fill.R
#
# A block is the text between
#   <!-- BEGIN generated: <name> -->  and  <!-- END generated -->
# and <name> is a path under dev/, optionally followed by
# " lines <from>-<to>" or " from <regex> to <regex>" to take part of the
# file. Markdown (a .md file, or one whose first line starts with "|" or
# "### ") is pasted as is; anything else goes in a text fence.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
f <- file.path(root, "remeasure-findings.md")
x <- readLines(f, warn = FALSE)
b <- grep("^<!-- BEGIN generated: .* -->$", x)
e <- grep("^<!-- END generated -->$", x)
stopifnot(length(b) == length(e), all(b < e))
out <- character()
prev <- 0L
for (k in seq_along(b)) {
  spec <- sub("^<!-- BEGIN generated: (.*) -->$", "\\1", x[b[k]])
  path <- sub(" (lines|from) .*$", "", spec)
  src <- readLines(file.path(root, path), warn = FALSE)
  if (grepl(" lines [0-9]+-[0-9]+$", spec)) {
    r <- as.integer(strsplit(sub("^.* lines ", "", spec), "-")[[1]])
    src <- src[r[1]:min(r[2], length(src))]
  } else if (grepl(" from .* to .*$", spec)) {
    fr <- sub("^.* from (.*) to (.*)$", "\\1", spec)
    to <- sub("^.* from (.*) to (.*)$", "\\2", spec)
    i <- grep(fr, src)[1]
    j <- grep(to, src); j <- j[j > i][1]
    stopifnot(!is.na(i), !is.na(j))
    src <- src[i:(j - 1)]
  }
  # a reducer's plain-text "## section" lines must not become headings
  # of the record, so only a .md file, a table or a "###" block is
  # pasted as markdown
  md <- grepl("[.]md$", path) ||
    (length(src) && grepl("^([|]|### )", src[1]))
  body <- if (md) src else c("```text", src, "```")
  out <- c(out, x[(prev + 1):b[k]], body)
  prev <- e[k] - 1L
}
out <- c(out, x[(prev + 1):length(x)])
writeLines(out, f, useBytes = TRUE)
cat("filled", length(b), "blocks in", f, "\n")

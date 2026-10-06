# Fill every generated block of the given records from
# dev/vigport-blocks.txt, so their counts are pasted and never typed.
#
#   Rscript dev/vigport-blocks.R > dev/vigport-blocks.txt
#   Rscript dev/vigport-fill.R <record.md> ...
#
# A record marks a block as
#   <!-- BEGIN generated: <name> -->
#   (anything)
#   <!-- END generated -->
# and the whole region is replaced by the block of that name. Run it
# again after any re-measurement; a name with no block is an error, so
# a renamed block cannot leave a stale copy behind.
root <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
src <- readLines(file.path(root, "vigport-blocks.txt"), warn = FALSE)
beg <- grep("^<!-- BEGIN generated: ", src)
blocks <- list()
for (a in beg) {
  nm <- sub("^<!-- BEGIN generated: (.*) -->$", "\\1", src[a])
  b <- a + grep("^<!-- END generated -->$", src[(a + 1):length(src)])[1]
  blocks[[nm]] <- src[a:b]
}
for (f in commandArgs(trailingOnly = TRUE)) {
  x <- readLines(f, warn = FALSE)
  out <- character()
  i <- 1L
  n <- 0L
  while (i <= length(x)) {
    if (grepl("^<!-- BEGIN generated: ", x[i])) {
      nm <- sub("^<!-- BEGIN generated: (.*) -->$", "\\1", x[i])
      if (is.null(blocks[[nm]])) stop(f, ": no generated block '", nm, "'")
      j <- i + grep("^<!-- END generated -->$", x[(i + 1):length(x)])[1]
      if (is.na(j)) stop(f, ": block '", nm, "' has no END line")
      out <- c(out, blocks[[nm]])
      i <- j + 1L
      n <- n + 1L
    } else {
      out <- c(out, x[i])
      i <- i + 1L
    }
  }
  writeLines(out, f)
  cat(f, ": filled", n, "blocks\n")
}

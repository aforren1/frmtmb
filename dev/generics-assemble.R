# Splice the generated measurement blocks into
# dev/generics-findings.md, replacing each <!--BLOCK:name--> marker
# with the captured output of the summarizer that produced it.
#
# The lane rule this implements: emit counts from the summarizer into a
# marked block and paste that block verbatim, so a figure cannot be
# typed and cannot go stale. A marker whose file is missing is replaced
# with a line that SAYS so, rather than left to look like prose.
root <- "C:/Users/adf44/source/r/frmtmb-wt-generics"
dir <- file.path(root, "dev", "generics-out2", "blocks")
md <- file.path(root, "dev", "generics-findings.md")
txt <- readLines(md, warn = FALSE)
marks <- grep("^<!--BLOCK:", txt)
out <- character()
for (i in seq_along(txt)) {
  if (!(i %in% marks)) { out <- c(out, txt[i]); next }
  nm <- sub("^<!--BLOCK:(.*)-->$", "\\1", txt[i])
  f <- file.path(dir, paste0(nm, ".txt"))
  if (!file.exists(f)) {
    out <- c(out, paste0("**BLOCK `", nm, "` NOT GENERATED.**"))
    cat("MISSING", nm, "\n")
  } else {
    out <- c(out, readLines(f, warn = FALSE))
    cat("spliced", nm, "\n")
  }
}
# writeBin, not writeLines: a text connection on Windows writes CRLF,
# which is how this file first came out CRLF in an LF tree
writeBin(charToRaw(paste0(paste(out, collapse = "\n"), "\n")), md)
cat("DONE\n")

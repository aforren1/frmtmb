# The one-substitution edit helper every dev/generics-edit*.R uses from
# round 2 on.
#
# It writes with writeBin, NOT writeLines.  On Windows `writeLines()`
# opens a text connection that turns every "\n" into "\r\n", so the
# round-1 edit scripts silently converted eight LF source files to
# CRLF in a tree that is LF.  See dev/generics-eol.R for the repair and
# for the worse thing the first repair did.
sub1 <- function(path, old, new, root =
                   "C:/Users/adf44/source/r/frmtmb-wt-generics") {
  f <- file.path(root, path)
  txt <- rawToChar(readBin(f, "raw", file.info(f)$size))
  if (!grepl(old, txt, fixed = TRUE)) stop("no match in ", path)
  if (length(gregexpr(old, txt, fixed = TRUE)[[1]]) != 1L) {
    stop("more than one match in ", path)
  }
  writeBin(charToRaw(sub(old, new, txt, fixed = TRUE)), f)
  cat("edited", path, "\n")
  invisible(TRUE)
}

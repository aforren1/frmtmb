# Realign the continuation lines of renamed function definitions.
for (f in commandArgs(trailingOnly = TRUE)) {
  x <- readLines(f)
  n <- 0L
  for (i in grep("^curve_[a-z_]+ <- function[(]", x)) {
    open <- regexpr("function[(]", x[i]) + 8L
    j <- i
    while (!grepl("[{] *$", x[j])) {
      j <- j + 1L
      x[j] <- paste0(strrep(" ", open), sub("^ *", "", x[j]))
      n <- n + 1L
    }
  }
  writeLines(x, f)
  cat(f, n, "\n")
}

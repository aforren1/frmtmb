# Restore the blank lines an earlier edit dropped from a test file: for each
# blank run in HEAD between lines A and B that now sit next to each other,
# put the blank line back. Usage: Rscript dev/impsamp-blank.R <path>
p <- commandArgs(TRUE)[[1]]
old <- system2("git", c("show", paste0("HEAD:", p)), stdout = TRUE)
cur <- readLines(p)
pairs <- list()
for (i in seq_along(old)) {
  if (nzchar(old[[i]]) || i == 1L || i == length(old)) next
  if (!nzchar(old[[i - 1L]]) || !nzchar(old[[i + 1L]])) next
  pairs[[length(pairs) + 1L]] <- c(old[[i - 1L]], old[[i + 1L]])
}
key <- vapply(pairs, paste, "", collapse = "\r")
out <- character(0)
n <- 0L
for (k in seq_along(cur)) {
  out <- c(out, cur[[k]])
  if (k < length(cur) && nzchar(cur[[k]]) && nzchar(cur[[k + 1L]]) &&
      paste(cur[[k]], cur[[k + 1L]], sep = "\r") %in% key) {
    out <- c(out, "")
    n <- n + 1L
  }
}
writeLines(out, p)
cat("restored", n, "blank lines\n")

# Sentences over 25 words in the top NEWS section of each package whose
# version changes at 0.72.0. Code spans count as one word. A sentence
# ends at ".", "!", "?" or ":" followed by a space, whatever the next
# character is: the release review's dev/relrev071-news.R splits only
# before a capital, so it also joins a sentence that starts with "brms"
# or a number to the one before it.
#   Rscript dev/rel072-news-long.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-release/"
files <- c(frmtmb = "NEWS.md", sample = "extensions/frmtmb.sample/NEWS.md",
           eam = "extensions/frmtmb.eam/NEWS.md",
           latent = "extensions/frmtmb.latent/NEWS.md",
           coupling = "extensions/frmtmb.coupling/NEWS.md")
for (k in names(files)) {
  L <- readLines(paste0(root, files[[k]]), warn = FALSE)
  h <- grep("^# ", L)
  sec <- L[h[1]:(h[2] - 1L)]
  starts <- grepl("^[*-] |^  - |^#|^$", sec)
  blocks <- split(sec, cumsum(starts))
  n <- 0
  for (bl in blocks) {
    t <- paste(trimws(bl), collapse = " ")
    if (!nzchar(t) || startsWith(t, "#")) next
    t <- gsub("`[^`]*`", "CODE", sub("^[*-] ", "", t))
    t <- gsub("**", "", t, fixed = TRUE)
    for (s in unlist(strsplit(t, "(?<=[.!?:])\\s+", perl = TRUE))) {
      w <- length(strsplit(trimws(s), "\\s+")[[1]])
      if (w > 25) {
        n <- n + 1
        cat(sprintf("  %s [%d] %s\n", k, w, s))
      }
    }
  }
  cat(k, ": sentences over 25 words:", n, "\n")
}

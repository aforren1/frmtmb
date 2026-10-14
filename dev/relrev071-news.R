# Release reviewer of 0.71.0: per package, the top NEWS section's
# headings, bullet count, and every sentence over 25 words (code spans
# count as one word each), plus em dashes and British spellings.
root <- "C:/Users/adf44/source/r/frmtmb-wt-release/"
files <- c(frmtmb = "NEWS.md",
           sample = "extensions/frmtmb.sample/NEWS.md",
           eam = "extensions/frmtmb.eam/NEWS.md",
           latent = "extensions/frmtmb.latent/NEWS.md",
           spline = "extensions/frmtmb.spline/NEWS.md")
for (k in names(files)) {
  L <- readLines(paste0(root, files[[k]]), warn = FALSE)
  h <- grep("^# ", L)
  sec <- L[h[1]:(h[2] - 1L)]
  cat("\n==", k, ":", sec[1], "| next:", L[h[2]], "\n")
  cat("  subheadings:", paste(grep("^## ", sec, value = TRUE),
                              collapse = " / "), "\n")
  b <- grep("^[*-] ", sec)
  cat("  bullets:", length(b), "\n")
  txt <- paste(sub("^[*-] ", "", sec[!grepl("^#", sec)]), collapse = " ")
  # split into bullets again to keep sentences inside one bullet
  starts <- c(b, length(sec) + 1L)
  long <- 0
  for (i in seq_along(b)) {
    bt <- paste(trimws(sec[b[i]:(starts[i + 1L] - 1L)]), collapse = " ")
    bt <- gsub("`[^`]*`", "CODE", sub("^[*-] ", "", bt))
    bt <- gsub("**", "", bt, fixed = TRUE)
    # sentences, and the items of a list inside a bullet
    s <- unlist(strsplit(bt, "(?<=[.!?:])\\s+(?=[A-Z(`])|\\s+-\\s+",
                         perl = TRUE))
    for (x in s) {
      n <- length(strsplit(trimws(x), "\\s+")[[1]])
      if (n > 25) {
        long <- long + 1
        cat(sprintf("  [%d words] %s\n", n, substr(x, 1, 160)))
      }
    }
  }
  cat("  sentences over 25 words:", long, "\n")
  cat("  em dashes:", sum(grepl("\u2014", sec)), "| BREAKING marks:",
      sum(grepl("BREAKING|[Bb]reaking", sec)), "\n")
}

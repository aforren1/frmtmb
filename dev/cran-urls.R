# Lane cran: urlchecker over the worktree (needs pandoc on PATH; reads
# public URLs only).
#
#   Rscript dev/cran-urls.R <worktree> > dev/cran-log/urls.txt
.libPaths(c("C:/Users/adf44/source/r/wt-cran-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- commandArgs(trailingOnly = TRUE)[1]
r <- urlchecker::url_check(root, progress = FALSE)
cat("urls with a problem:", nrow(r), "\n")
if (nrow(r)) {
  df <- as.data.frame(r)
  for (j in seq_len(nrow(df))) {
    cat("URL\t", df$URL[j], "\t", df$Status[j], "\t", df$Message[j], "\t",
        df$New[j], "\t", paste(unlist(df$From[j]), collapse = ","), "\n",
        sep = "")
  }
}

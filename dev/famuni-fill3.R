# Replace fill2's FINAL block function with the complete provenance list.
p <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-fill2.R"
x <- paste(readLines(p), collapse = "\n")
old <- '  fs <- c("final-installs.txt", "final-suite.txt", "final-gated.txt",
          "final-drmtmb.txt", "final-openblas.txt")'
new <- '  fs <- c("final-installs.txt", "final-install-sample2.txt",
          "final-suite.txt", "final-gated.txt", "final-drmtmb.txt",
          "final-reach.txt", "final-openblas.txt")'
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)

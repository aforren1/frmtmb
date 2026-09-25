# List each replicate's error or diagnose() verdict, for one arm.
arm <- commandArgs(trailingOnly = TRUE)[[1]]
for (f in list.files("dev/phase3b-log/recov", paste0("^", arm, "-"),
                     full.names = TRUE)) {
  x <- readRDS(f)
  cat(basename(f), ":", if (!is.null(x$error)) paste("ERROR", x$error) else
    paste(x$diagnose, collapse = " | "), "\n")
}

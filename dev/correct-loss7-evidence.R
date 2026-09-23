# The SEVENTH loss of the user library, 2026-09-22. Evidence collected
# BEFORE the restore, in the order dev/machine-library.md asks for it:
# the full sorted mtime list, whether anything outside %LOCALAPPDATA%
# was hit, whether a 00LOCK is present, and the canary.
#   Rscript dev/correct-loss7-evidence.R
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
others <- c("C:/Users/adf44/source/r/correct-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib")
hollow_of <- function(lib) {
  d <- list.dirs(lib, recursive = FALSE, full.names = FALSE)
  d <- setdiff(d, grep("^00LOCK", d, value = TRUE))
  list(dirs = d, hollow = d[!file.exists(file.path(lib, d, "DESCRIPTION"))])
}
u <- hollow_of(ul)
cat("user library:", length(u$dirs), "dirs,", length(u$hollow), "hollow\n")
mt <- file.info(file.path(ul, u$hollow))$mtime
o <- order(mt)
cat("mtime range:", format(min(mt)), "to", format(max(mt)), "\n")
cat("the full sorted list, one per line:\n")
for (i in o) cat(sprintf("  %s  %s\n", format(mt[i], "%Y-%m-%d %H:%M:%S"),
                         u$hollow[i]))
cat("\nseconds between first and last:",
    as.numeric(difftime(max(mt), min(mt), units = "secs")), "\n")
cat("00LOCK directories present:",
    paste(grep("^00LOCK", list.dirs(ul, recursive = FALSE,
                                    full.names = FALSE), value = TRUE),
          collapse = " "), "\n")
cat("canary:", file.exists(file.path(ul, "ZZZ-canary.txt")), "\n")
if (file.exists(file.path(ul, "ZZZ-canary.txt"))) {
  cat("canary mtime:",
      format(file.info(file.path(ul, "ZZZ-canary.txt"))$mtime), "\n")
}
for (lib in others) {
  h <- hollow_of(lib)
  cat(sprintf("outside %%LOCALAPPDATA%%: %-45s dirs %d hollow %d\n", lib,
              length(h$dirs), length(h$hollow)))
}
cat("\nsurvivors, for the record:",
    paste(setdiff(u$dirs, u$hollow), collapse = " "), "\n")

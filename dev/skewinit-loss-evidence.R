# Evidence for the seventh library loss, collected BEFORE any restore,
# in the order dev/machine-library.md asks for it.
libs <- c(user = "C:/Users/adf44/AppData/Local/R/win-library/4.6",
          rellib = "C:/Users/adf44/source/r/rellib-r3",
          pinlib = "C:/Users/adf44/source/r/pinlib",
          lane = "C:/Users/adf44/source/r/skewinit-lib")
for (nm in names(libs)) {
  L <- libs[[nm]]
  if (!dir.exists(L)) { cat(nm, ": MISSING\n", sep = ""); next }
  d <- list.dirs(L, recursive = FALSE, full.names = FALSE)
  has <- file.exists(file.path(L, d, "DESCRIPTION"))
  cat(sprintf("%-7s %s\n", nm, L))
  cat(sprintf("        %d dirs, %d with DESCRIPTION, %d hollow\n",
              length(d), sum(has), sum(!has)))
  lock <- grep("^00LOCK", d, value = TRUE)
  cat("        00LOCK dirs:", if (length(lock)) paste(lock, collapse = " ")
      else "none", "\n")
  if (any(!has)) {
    mt <- file.info(file.path(L, d[!has]))$mtime
    o <- order(mt)
    cat("        hollow mtimes, sorted:\n")
    print(table(format(sort(mt), "%Y-%m-%d %H:%M")))
    cat("        first:", format(min(mt)), " last:", format(max(mt)), "\n")
    cat("        first five hollow:",
        paste(utils::head(d[!has][o], 5), collapse = " "), "\n")
  }
  cy <- file.path(L, "ZZZ-canary.txt")
  if (file.exists(cy)) {
    cat("        ZZZ-canary.txt SURVIVED, mtime",
        format(file.info(cy)$mtime), "\n")
  } else if (nm == "user") {
    cat("        ZZZ-canary.txt DESTROYED\n")
  }
}

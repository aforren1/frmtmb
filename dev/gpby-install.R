# Roxygenise and install the lane's core (and named extensions) into the
# lane library. Usage: Rscript dev/gpby-install.R [core] [ext ...]
LIB <- "C:/Users/adf44/source/r/wt-gpby-lib"
WT <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
args <- commandArgs(TRUE)
if (!length(args)) args <- "core"
rcmd <- file.path(R.home("bin"), "R.exe")
for (a in args) {
  path <- if (a == "core") WT else file.path(WT, "extensions", a)
  cat("== roxygenise", path, "\n")
  roxygen2::roxygenise(path)
  cat("== install", path, "\n")
  st <- system2(rcmd, c("CMD", "INSTALL", paste0("--library=", LIB),
                        "--no-multiarch", "--no-test-load", shQuote(path)),
                stdout = TRUE, stderr = TRUE)
  cat(tail(st, 5), sep = "\n")
  if (!is.null(attr(st, "status")) && attr(st, "status") != 0) {
    cat(st, sep = "\n")
    stop("install failed: ", a)
  }
}
cat("INSTALLED\n")

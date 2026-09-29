# Reviewer, claim 7: knit vignette("case-studies") on the worker's build,
# every chunk, with knitr (no pandoc needed), and report the first error.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
cat("frmtmb:", find.package("frmtmb"), "\n")
src <- "C:/Users/adf44/source/r/frmtmb-wt-resmooth/vignettes/case-studies.Rmd"
wd <- file.path(tempdir(), "revvig")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)
file.copy(src, file.path(wd, "case-studies.Rmd"), overwrite = TRUE)
old <- setwd(wd)
knitr::opts_chunk$set(error = FALSE)
t0 <- proc.time()
r <- tryCatch(knitr::knit("case-studies.Rmd", quiet = TRUE),
              error = function(e) e)
el <- (proc.time() - t0)[["elapsed"]]
setwd(old)
if (inherits(r, "error")) {
  cat("KNIT FAILED after", sprintf("%.0f", el), "s:\n")
  cat(conditionMessage(r), "\n")
} else {
  md <- readLines(file.path(wd, r), warn = FALSE)
  cat("KNIT OK in", sprintf("%.0f", el), "s | output lines:", length(md),
      "\n")
  cat("figures written:",
      length(list.files(wd, pattern = "\\.(png|svg)$", recursive = TRUE)),
      "\n")
  cat("lines matching 'Error':", length(grep("^## Error", md)), "\n")
  if (length(grep("^## Error", md))) {
    cat(paste(grep("^## Error", md, value = TRUE), collapse = "\n"), "\n")
  }
}

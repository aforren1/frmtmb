# one test file per process; usage: Rscript es-runtest.R <file>
lib_es <- "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/529b6e73-d28f-46aa-a279-7dbeeb58fd4f/scratchpad/es-lib"
.libPaths(c(lib_es, "C:/Users/adf44/AppData/Local/R/win-library/4.6", .libPaths()))
library(testthat); library(frmtmb)
f <- commandArgs(trailingOnly = TRUE)[[1]]
wt <- "C:/Users/adf44/source/r/frmtmb-wt-esicar"
setwd(file.path(wt, "tests", "testthat"))
r <- as.data.frame(test_file(f, reporter = "silent", package = "frmtmb"))
cat(sprintf("FILE %s tests %d pass %d fail %d error %d skip %d\n",
            f, nrow(r), sum(r$passed), sum(r$failed),
            sum(r$error), sum(r$skipped)))
bad <- r[r$failed > 0 | r$error, , drop = FALSE]
if (nrow(bad)) {
  for (i in seq_len(nrow(bad))) cat("  BAD:", bad$test[i], "\n")
  print(utils::head(as.data.frame(r)[r$failed > 0 | r$error, c("test", "failed", "error")]))
}
